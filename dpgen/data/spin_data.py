"""Stage 5: collect the final frame of each magnetic VASP task.

DPData owns the ordinary DeepMD labels.  Magnetic labels are read from the
last OUTCAR/OSZICAR blocks and written alongside DPData's raw/npy output.
"""

import json
import re
import tempfile
from pathlib import Path

import dpdata
import numpy as np

from dpgen import dlog
from dpgen.data.spin_tasks import _saved_tasks, check_spin_results, read_spin_incar

DATA_DIR = "04.data"
RMSE_LIMIT = 5.0e-3
_ATOM_ROW = re.compile(r"^\s*(\d+)\s+(.*)$")
_MAGNETIZATION = re.compile(r"^\s*magnetization\s*\(([xyz])\)\s*$", re.I)
_VASP_VERSION = re.compile(r"\bvasp\.(\d+)(?:\.\d+)*\b", re.I)


def _number(token, context):
    try:
        value = float(token.replace("D", "E").replace("d", "e"))
    except ValueError as error:
        raise ValueError(f"{context}: invalid number {token!r}") from error
    if not np.isfinite(value):
        raise ValueError(f"{context}: nonfinite number {token!r}")
    return value


def _outcar_magnetization_block(lines, marker, natoms, path):
    """Read one OUTCAR orbital table, taking its final ``tot`` column."""
    values = {}
    header_seen = False
    for line in lines[marker + 1 :]:
        if _MAGNETIZATION.match(line):
            break
        if "# of ion" in line.lower():
            header_seen = True
            if "tot" not in line.lower().split():
                raise ValueError(f"{path}: magnetization table has no tot column")
            continue
        if not header_seen:
            continue
        match = _ATOM_ROW.match(line)
        if match is None:
            if values:
                break
            continue
        index = int(match.group(1))
        tokens = match.group(2).split()
        if index < 1 or index > natoms or index in values or not tokens:
            raise ValueError(f"{path}: invalid magnetization atom row {line.strip()!r}")
        values[index] = _number(tokens[-1], str(path))
        if len(values) == natoms:
            break
    if sorted(values) != list(range(1, natoms + 1)):
        raise ValueError(
            f"{path}: incomplete magnetization table; found atoms "
            f"{sorted(values)}, expected 1..{natoms}"
        )
    return np.array([values[index] for index in range(1, natoms + 1)])


def read_final_magnetization(outcar, natoms):
    """Return VASP major version and the final noncollinear moments (N, 3)."""
    outcar = Path(outcar)
    lines = outcar.read_text(errors="replace").splitlines()
    version = next(
        (match for line in lines if (match := _VASP_VERSION.search(line))), None
    )
    if version is None or int(version.group(1)) not in (5, 6):
        raise ValueError(f"{outcar}: cannot determine supported VASP 5/6 version")
    markers = [
        (index, match.group(1).lower())
        for index, line in enumerate(lines)
        if (match := _MAGNETIZATION.match(line))
    ]
    last_x = next((index for index, axis in reversed(markers) if axis == "x"), None)
    if last_x is None:
        raise ValueError(f"{outcar}: missing magnetization (x) block")
    final = [(index, axis) for index, axis in markers if index >= last_x]
    if [axis for _, axis in final] != ["x", "y", "z"]:
        raise ValueError(
            f"{outcar}: final noncollinear magnetization requires x/y/z blocks"
        )
    columns = [
        _outcar_magnetization_block(lines, index, natoms, outcar)
        for index, _ in final
    ]
    return int(version.group(1)), np.column_stack(columns)


def _oszicar_rows(lines, marker, natoms, path):
    """Read indexed vector rows; unlisted unconstrained atoms remain zero."""
    values = np.zeros((natoms, 3), dtype=float)
    found = set()
    for line in lines[marker + 1 :]:
        if re.search(r"\bMW_int\b|lambda\s*\*\s*MW_perp", line):
            break
        match = _ATOM_ROW.match(line)
        if match is None:
            if found:
                break
            continue
        index = int(match.group(1))
        tokens = match.group(2).split()
        if len(tokens) < 3 or index < 1 or index > natoms or index in found:
            if found:
                break
            raise ValueError(f"{path}: invalid magnetic vector row {line.strip()!r}")
        values[index - 1] = [_number(token, str(path)) for token in tokens[:3]]
        found.add(index)
        if len(found) == natoms:
            break
    if not found:
        raise ValueError(f"{path}: magnetic vector section has no atom rows")
    return values


def read_final_spin_force(oszicar, natoms, vasp_major):
    """Return magnetic force in eV, using the VASP 5/6 factor convention."""
    oszicar = Path(oszicar)
    if vasp_major not in (5, 6):
        raise ValueError(f"{oszicar}: unsupported VASP major version {vasp_major}")
    lines = oszicar.read_text(errors="replace").splitlines()
    mw_markers = [i for i, line in enumerate(lines) if re.search(r"\bMW_int\b", line)]
    force_markers = [
        i for i, line in enumerate(lines) if re.search(r"lambda\s*\*\s*MW_perp", line)
    ]
    if not mw_markers or not force_markers:
        raise ValueError(
            f"{oszicar}: missing MW_int or lambda*MW_perp magnetic section"
        )
    mw_index, force_index = mw_markers[-1], force_markers[-1]
    if mw_index >= force_index:
        raise ValueError(
            f"{oszicar}: final MW_int/lambda*MW_perp sections out of order"
        )
    mw = _oszicar_rows(lines, mw_index, natoms, oszicar)
    perpendicular = _oszicar_rows(lines, force_index, natoms, oszicar)
    factor = 2.0 if vasp_major == 5 else 1.0
    return factor * perpendicular * np.linalg.norm(mw, axis=1)[:, None]


def magnetic_rmse(initial, final):
    """RMSE of moment magnitudes on atoms with nonzero initial moments."""
    initial_length = np.linalg.norm(initial, axis=1)
    mask = initial_length > 0
    if not np.any(mask):
        raise ValueError("INCAR has no nonzero initial magnetic moments")
    final_length = np.linalg.norm(final, axis=1)
    return float(np.sqrt(np.mean((initial_length[mask] - final_length[mask]) ** 2)))


def _last_dpdata_frame(outcar):
    """Read the last ionic step, including a normally terminated unconverged step.

    DPData rotates VASP cells, coordinates, and atomic forces to a lower-
    triangular cell.  Read the unrotated final cell too, so magnetic vectors
    can be transformed by the same orthogonal matrix.
    """
    system = dpdata.LabeledSystem(
        str(outcar), fmt="vasp/outcar", convergence_check=False
    )
    if system.get_nframes() == 0 or not system.has_forces():
        raise ValueError(f"{outcar}: DPData found no labeled ionic frames")
    force_blocks = outcar.read_text(errors="replace").count("TOTAL-FORCE")
    if force_blocks != system.get_nframes():
        raise ValueError(
            f"{outcar}: {force_blocks} VASP force blocks but DPData read "
            f"{system.get_nframes()} frames; final-frame alignment is uncertain"
        )
    frame = system.sub_system(-1)
    raw = dpdata.vasp.outcar.get_frames(str(outcar), convergence_check=False)
    raw_cells = raw[3]
    raw_coords = raw[4]
    if len(raw_cells) != system.get_nframes() or len(raw_coords) != len(raw_cells):
        raise ValueError(f"{outcar}: DPData raw/rotated frame counts differ")
    raw_cell = np.asarray(raw_cells[-1], dtype=float)
    rotated_cell = np.asarray(frame["cells"][0], dtype=float)
    rotation = np.linalg.solve(raw_cell, rotated_cell)
    if not np.allclose(rotation.T @ rotation, np.eye(3), atol=1e-6) or not np.isclose(
        np.linalg.det(rotation), 1.0, atol=1e-6
    ):
        raise ValueError(f"{outcar}: cannot align magnetic and DPData coordinates")
    if not np.allclose(
        np.asarray(raw_coords[-1], dtype=float) @ rotation,
        frame["coords"][0],
        atol=1e-5,
    ):
        raise ValueError(f"{outcar}: final DPData coordinates do not match VASP frame")
    for name in ("cells", "coords", "energies", "forces"):
        if name not in frame.data or not np.isfinite(frame[name]).all():
            raise ValueError(f"{outcar}: missing or nonfinite final {name}")
    return frame, rotation, system.get_nframes() - 1


def _group_name(frame, groups):
    """Separate differing atom orders instead of letting DPData reorder spins."""
    formula = frame.formula
    atom_order = tuple(int(value) for value in frame["atom_types"])
    key = (formula, atom_order)
    if key not in groups:
        used = {group["name"] for group in groups.values()}
        name = formula
        suffix = 1
        while name in used:
            name = f"{formula}-layout-{suffix:03d}"
            suffix += 1
        groups[key] = {"name": name, "system": frame.copy(), "records": []}
    else:
        groups[key]["system"].append(frame)
    return groups[key]


def collect_spin_data(jdata):
    """Filter Stage 4 tasks and atomically create aligned DeepMD raw/npy data."""
    stage, tasks = _saved_tasks(jdata)
    output_root = Path(jdata.get("out_dir", ".")).resolve()
    destination = output_root / DATA_DIR
    if destination.exists() or destination.is_symlink():
        raise RuntimeError(f"output stage already exists: {destination}")
    # An incomplete VASP task is not an RMSE rejection.  Keep the existing
    # Stage 4 completion check separate from magnetic data parsing/filtering.
    check_spin_results(jdata)

    groups = {}
    rejected = []
    for task in tasks:
        task_dir = stage / task
        outcar = task_dir / "OUTCAR"
        oszicar = task_dir / "OSZICAR"
        try:
            from pymatgen.io.vasp.inputs import Poscar

            frame, rotation, frame_index = _last_dpdata_frame(outcar)
            natoms = frame.get_natoms()
            poscar_species = [
                site.specie.symbol
                for site in Poscar.from_file(task_dir / "POSCAR").structure
            ]
            dpdata_species = [
                frame["atom_names"][int(atom_type)]
                for atom_type in frame["atom_types"]
            ]
            if poscar_species != dpdata_species:
                raise ValueError(
                    f"{outcar}: DPData atom order/species differs from "
                    f"{task_dir / 'POSCAR'}"
                )
            _, initial = read_spin_incar(task_dir / "INCAR", natoms)
            major, final = read_final_magnetization(outcar, natoms)
            magnetic_force = read_final_spin_force(oszicar, natoms, major)
            rmse = magnetic_rmse(initial, final)
        except Exception as error:
            raise RuntimeError(f"spin data task {task}: {error}") from error
        if rmse > RMSE_LIMIT:
            rejected.append(
                {"task": task, "rmse": rmse, "reason": "RMSE exceeds 5.0e-3"}
            )
            dlog.warning(
                "spin data task %s rejected: RMSE %.8g > %.8g",
                task,
                rmse,
                RMSE_LIMIT,
            )
            continue

        # The standard DeepMD data and these magnetic vectors must use the
        # same Cartesian basis and the same atom order.
        lengths = np.linalg.norm(final, axis=1)
        directions = np.divide(
            final,
            lengths[:, None],
            out=np.zeros_like(final),
            where=lengths[:, None] != 0,
        ) @ rotation
        magnetic_force = magnetic_force @ rotation
        frame.sort_atom_names()
        group = _group_name(frame, groups)
        group["records"].append(
            {
                "task": task,
                "frame": frame_index,
                "rmse": rmse,
                "spin": directions,
                "spin_force": magnetic_force,
                "spin_length": lengths,
            }
        )

    if not groups:
        raise RuntimeError(
            f"no spin tasks passed RMSE <= {RMSE_LIMIT:g}; "
            f"rejected {len(rejected)} of {len(tasks)} tasks"
        )

    with tempfile.TemporaryDirectory(dir=str(output_root)) as temporary:
        scratch = Path(temporary) / DATA_DIR
        deepmd = scratch / "deepmd"
        deepmd.mkdir(parents=True)
        selected = []
        for group in groups.values():
            folder = deepmd / group["name"]
            system = group["system"]
            records = group["records"]
            count = len(records)
            system.to_deepmd_raw(str(folder))
            system.to_deepmd_npy(str(folder), set_size=count)
            set_folder = folder / "set.000"
            for field in ("spin", "spin_force", "spin_length"):
                values = np.stack([record[field] for record in records]).reshape(
                    count, -1
                )
                np.savetxt(folder / f"{field}.raw", values)
                np.save(set_folder / f"{field}.npy", values.astype(np.float32))
            frames = [
                {
                    "index": index,
                    "task": record["task"],
                    "outcar_frame": record["frame"],
                    "rmse": record["rmse"],
                }
                for index, record in enumerate(records)
            ]
            (folder / "frames.json").write_text(json.dumps(frames, indent=2))
            selected.extend({"system": group["name"], **record} for record in frames)
        (deepmd / "selection.json").write_text(
            json.dumps(
                {
                    "rmse_limit": RMSE_LIMIT,
                    "selected": selected,
                    "rejected": rejected,
                },
                indent=2,
            )
        )
        scratch.replace(destination)
    return len(selected)
