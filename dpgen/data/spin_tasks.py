"""Stage 4: noncollinear VASP tasks for each exported AIMD snapshot.

The baseline is always generated; a perturbation provider supplies additional
named magnetic configurations such as C1 and C2.
"""

import filecmp
import json
import re
import tempfile
from pathlib import Path

import numpy as np

from dpgen.data.spin_init import (
    _expected_task_paths,
    _make_relative_symlink,
    _require_file,
    _unique,
)
from dpgen.dispatcher.Dispatcher import make_submission
from dpgen.generator.lib.utils import check_api_version

SPIN_DIR = "03.spin"


def _vectors(values, natoms, context):
    try:
        result = np.asarray(values, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{context}: magnetic moments must be numeric") from error
    if result.size != 3 * natoms or not np.isfinite(result).all():
        raise ValueError(f"{context}: expected {3 * natoms} finite magnetic components")
    if result.ndim > 2 or (result.ndim == 2 and result.shape != (natoms, 3)):
        raise ValueError(f"{context}: expected shape ({natoms}, 3)")
    return result.reshape(natoms, 3).copy()


def _read_vector_tag(value, natoms, context):
    # pymatgen versions used with Python 3.9 leave M_CONSTR as a string and
    # some versions parse MAGMOM scientific notation incorrectly. Only these
    # two numeric fields need this adapter; other INCAR tags use pymatgen.
    result = []
    try:
        for token in value.split():
            count, number = token.split("*", 1) if "*" in token else ("1", token)
            if not count.isdigit() or int(count) <= 0:
                raise ValueError("invalid repetition")
            if len(result) + int(count) > 3 * natoms:
                raise ValueError("too many components")
            result.extend(
                [float(number.replace("D", "E").replace("d", "e"))] * int(count)
            )
    except ValueError as error:
        raise ValueError(f"{context}: invalid vector: {error}") from error
    return _vectors(result, natoms, context)


def read_spin_incar(path, natoms):
    """Read a static noncollinear template with matching MAGMOM/M_CONSTR."""
    from pymatgen.io.vasp.inputs import Incar

    path = _require_file(path, "spin INCAR source")
    text = path.read_text()
    lines = [
        re.split(r"[!#]", line, maxsplit=1)[0].rstrip() for line in text.splitlines()
    ]
    if lines and lines[-1].endswith("\\"):
        raise ValueError(f"{path}: unfinished INCAR continuation")
    logical = "\n".join(lines).replace("\\\n", " ")
    fields = {}
    for assignment in re.split(r"[;\n]", logical):
        match = re.match(r"\s*(MAGMOM|M_CONSTR)\s*=\s*(.*)$", assignment, re.I)
        if match:
            key = match[1].upper()
            if key in fields:
                raise ValueError(f"{path}: duplicate {key}")
            fields[key] = _read_vector_tag(match[2], natoms, f"{path}: {key}")
    if set(fields) != {"MAGMOM", "M_CONSTR"}:
        raise ValueError(f"{path}: both MAGMOM and M_CONSTR are required")
    if not np.allclose(fields["MAGMOM"], fields["M_CONSTR"], rtol=0, atol=1e-12):
        raise ValueError(f"{path}: MAGMOM and M_CONSTR must match")
    try:
        incar = Incar.from_str(text)
        if not incar.get("LNONCOLLINEAR", incar.get("LSORBIT", False)):
            raise ValueError("three-component moments require LNONCOLLINEAR = .TRUE.")
        if int(incar.get("NSW", 0)) != 0 or int(incar.get("IBRION", -1)) != -1:
            raise ValueError("Stage 4 requires static calculations (NSW=0, IBRION=-1)")
    except Exception as error:
        raise ValueError(f"{path}: invalid spin INCAR: {error}") from error
    return incar, fields["MAGMOM"]


def write_spin_incar(template, moments, path):
    """Write a separate INCAR with identical initial and target vectors."""
    values = _vectors(moments, len(moments), str(path)).reshape(-1)
    incar = template.copy()
    incar.pop("MAGMOM", None)
    incar.pop("M_CONSTR", None)
    # Avoid version-specific MAGMOM compression/re-parsing when writing vectors.
    vector_text = " ".join(format(value, ".16g") for value in values)
    body = str(incar)
    if body and not body.endswith("\n"):
        body += "\n"
    Path(path).write_text(body + f"MAGMOM = {vector_text}\nM_CONSTR = {vector_text}\n")


def _forward_sources(mdata):
    sources = {}
    for value in mdata.get("fp_user_forward_files", []):
        name = Path(value).name
        if (
            name in {"POSCAR", "POTCAR", "INCAR", "OUTCAR", "OSZICAR", "fp.log"}
            or name in sources
        ):
            raise ValueError(f"reserved or duplicate spin forward filename: {name}")
        sources[name] = _require_file(value, "spin forward file")
    return sources


def _snapshot_poscars(jdata):
    """Discover every valid snapshot while requiring all configured parents."""
    source = Path(jdata.get("out_dir", ".")).resolve() / "02.disp"
    if not source.is_dir():
        raise FileNotFoundError(f"snapshot stage does not exist: {source}")

    discovered = {}
    for scale in sorted(
        (path for path in source.iterdir() if path.is_dir()),
        key=lambda path: path.name,
    ):
        if not re.fullmatch(r"scale-[0-9]+(?:\.[0-9]+)?", scale.name):
            if any(scale.rglob("POSCAR")):
                raise ValueError(f"invalid scale snapshot directory: {scale}")
            continue
        for perturbation in scale.iterdir():
            if not perturbation.is_dir():
                continue
            if not re.fullmatch(r"\d{6}", perturbation.name):
                if any(perturbation.rglob("POSCAR")):
                    raise ValueError(
                        f"invalid perturbation snapshot directory: {perturbation}"
                    )
                continue
            parent = f"{scale.name}/{perturbation.name}"
            frames = []
            for frame in perturbation.iterdir():
                if not frame.is_dir():
                    continue
                if not re.fullmatch(r"\d{2,}", frame.name):
                    if (frame / "POSCAR").exists():
                        raise ValueError(f"invalid frame snapshot directory: {frame}")
                    continue
                frames.append(
                    (
                        int(frame.name),
                        frame.name,
                        _require_file(
                            frame / "POSCAR",
                            "snapshot POSCAR",
                            f"{parent}/{frame.name}",
                        ),
                    )
                )
            if frames:
                discovered[parent] = [item[2] for item in sorted(frames)]

    missing = [
        parent for parent in _expected_task_paths(jdata) if parent not in discovered
    ]
    if missing:
        raise FileNotFoundError(
            f"no snapshots for configured task(s) {', '.join(missing)} under: {source}"
        )
    if not discovered:
        raise ValueError(f"no snapshots under {source}")

    def parent_key(value):
        scale, perturbation = value.split("/")
        return float(scale[6:]), int(perturbation)

    return [
        (parent, poscar)
        for parent in sorted(discovered, key=parent_key)
        for poscar in discovered[parent]
    ]


def plan_spin_tasks(jdata, mdata, perturb=None):
    """Validate all inputs and return one task per snapshot/configuration.

    ``perturb(moments, count)`` returns a mapping of safe task names to N x 3
    arrays. The baseline 000000 is always added here, not by the provider.
    """
    from pymatgen.io.vasp.inputs import Poscar

    count = jdata.get("spin_pert_numb", 0)
    if isinstance(count, bool) or not isinstance(count, int) or count < 0:
        raise ValueError("spin_pert_numb must be a non-negative integer")
    if count and perturb is None:
        raise NotImplementedError(
            "magnetic perturbation algorithm not supplied; use spin_pert_numb=0"
        )
    _forward_sources(mdata)
    if not jdata.get("potcars"):
        raise ValueError("potcars must contain at least one file")
    for path in jdata["potcars"]:
        _require_file(path, "spin POTCAR source")
    tasks = []
    species_reference = None
    for parent, poscar in _snapshot_poscars(jdata):
        frame_name = poscar.parent.name
        try:
            structure = Poscar.from_file(poscar).structure
        except Exception as error:
            raise ValueError(f"invalid snapshot {poscar}: {error}") from error
        species = [str(site.specie) for site in structure]
        if species_reference is not None and species != species_reference:
            raise ValueError(f"snapshot atom order/species differs: {poscar}")
        species_reference = species
        try:
            incar, moments = read_spin_incar(jdata["spin_incar"], len(structure))
        except FileNotFoundError as error:
            raise FileNotFoundError(f"{parent}/{frame_name}: {error}") from error
        except ValueError as error:
            raise ValueError(f"{parent}/{frame_name}: {error}") from error
        configs = {"000000": moments}
        if count:
            try:
                extra = perturb(moments.copy(), count)
            except ValueError as error:
                raise ValueError(f"{poscar}: {error}") from error
            if not isinstance(extra, dict) or len(extra) != count:
                raise ValueError(
                    f"{poscar}: perturbation provider must return {count} named arrays"
                )
            for name, values in extra.items():
                if (
                    not isinstance(name, str)
                    or not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_-]*", name)
                    or name == "000000"
                ):
                    raise ValueError(f"invalid magnetic configuration name: {name!r}")
                configs[name] = _vectors(values, len(structure), f"{poscar}: {name}")
        for name, values in configs.items():
            tasks.append(
                {
                    "task": f"{parent}/{frame_name}/{name}",
                    "poscar": poscar,
                    "incar": incar,
                    "moments": values,
                }
            )
    return tasks


def make_spin_tasks(jdata, mdata, perturb=None):
    """Atomically create 03.spin with real relative links and private INCARs."""
    root = Path(jdata.get("out_dir", ".")).resolve()
    stage = root / SPIN_DIR
    if stage.exists() or stage.is_symlink():
        raise RuntimeError(
            f"output stage already exists: {stage}; use spin_action=run to submit existing tasks"
        )
    tasks = plan_spin_tasks(jdata, mdata, perturb)
    forwards = _forward_sources(mdata)
    with tempfile.TemporaryDirectory(dir=str(root)) as temporary:
        scratch = Path(temporary) / SPIN_DIR
        scratch.mkdir()
        with (scratch / "POTCAR").open("wb") as output:
            for path in jdata["potcars"]:
                output.write(Path(path).read_bytes())
        for task in tasks:
            dest = scratch / task["task"]
            final = stage / task["task"]
            dest.mkdir(parents=True)
            _make_relative_symlink(dest / "POSCAR", task["poscar"], final / "POSCAR")
            _make_relative_symlink(dest / "POTCAR", stage / "POTCAR", final / "POTCAR")
            write_spin_incar(task["incar"], task["moments"], dest / "INCAR")
            for name, path in forwards.items():
                _make_relative_symlink(dest / name, path, final / name)
        paths = [task["task"] for task in tasks]
        (scratch / "tasks.json").write_text(json.dumps({"tasks": paths}, indent=2))
        scratch.replace(stage)
    return paths


def _saved_tasks(jdata):
    stage = Path(jdata.get("out_dir", ".")).resolve() / SPIN_DIR
    manifest = _require_file(stage / "tasks.json", "spin task manifest")
    try:
        data = json.loads(manifest.read_text())
        paths = data["tasks"]
    except (json.JSONDecodeError, KeyError, TypeError) as error:
        raise ValueError(f"{manifest}: invalid spin task manifest: {error}") from error
    if not isinstance(paths, list) or not paths or len(set(paths)) != len(paths):
        raise ValueError(f"{manifest}: expected unique nonempty task list")
    for task in paths:
        if not isinstance(task, str) or not re.fullmatch(
            r"scale-[0-9]+(?:\.[0-9]+)?/\d{6}/\d{2,}/[A-Za-z0-9][A-Za-z0-9_-]*",
            task,
        ):
            raise ValueError(f"{manifest}: invalid task path {task!r}")
        if stage not in (stage / task).resolve().parents:
            raise ValueError(f"{manifest}: task escapes stage: {task}")
    return stage, paths


def materialize_spin_forward_files(stage, paths, forwards):
    """Add missing user-forward links to an existing, validated task tree."""
    stage = Path(stage).resolve()
    pending = []
    for task in paths:
        for name, source in forwards.items():
            link = stage / task / name
            if link.exists() or link.is_symlink():
                existing = _require_file(link, "spin task forward file", task)
                if link.is_symlink():
                    matches = existing == Path(source).resolve()
                else:
                    matches = filecmp.cmp(existing, source, shallow=False)
                if not matches:
                    raise RuntimeError(
                        f"spin task {task} already has a different {name}: {link}"
                    )
            else:
                pending.append((link, source))

    created = []
    try:
        for link, source in pending:
            _make_relative_symlink(link, source)
            created.append(link)
    except Exception:
        for link in created:
            link.unlink(missing_ok=True)
        raise


def make_spin_submission(jdata, mdata):
    """Build real dpdispatcher objects; no job is submitted by this function."""
    from pymatgen.io.vasp.inputs import Poscar

    stage, paths = _saved_tasks(jdata)
    forwards = _forward_sources(mdata)
    for task in paths:
        dest = stage / task
        for name in ["POSCAR", "INCAR", "POTCAR"]:
            _require_file(dest / name, "spin task input", task)
        try:
            natoms = len(Poscar.from_file(dest / "POSCAR").structure)
        except Exception as error:
            raise ValueError(
                f"invalid POSCAR for spin task {task}: {dest / 'POSCAR'}: {error}"
            ) from error
        read_spin_incar(dest / "INCAR", natoms)
    materialize_spin_forward_files(stage, paths, forwards)
    check_api_version(mdata)
    context = mdata["fp_machine"].get("context_type", "").lower()
    if context in ("lebesguecontext", "dpcloudservercontext", "openapicontext"):
        try:
            import oss2  # noqa: F401
        except ImportError as error:
            raise RuntimeError(
                'Bohrium dependency unavailable; run: python -m pip install "dpdispatcher[bohrium]"'
            ) from error
    return make_submission(
        mdata["fp_machine"],
        mdata["fp_resources"],
        commands=[mdata["fp_command"]],
        work_path=str(stage),
        run_tasks=paths,
        group_size=mdata["fp_group_size"],
        forward_common_files=[],
        forward_files=["POSCAR", "INCAR", "POTCAR", *forwards],
        backward_files=_unique(
            ["OUTCAR", "OSZICAR"] + mdata.get("fp_user_backward_files", [])
        ),
        outlog="fp.log",
        errlog="fp.log",
    )


def check_spin_results(jdata):
    """Check static VASP termination/files, not magnetic RMSE or convergence."""
    stage, paths = _saved_tasks(jdata)
    for task in paths:
        outcar = _require_file(stage / task / "OUTCAR", "spin OUTCAR", task)
        text = outcar.read_text(errors="replace")
        if not text:
            raise RuntimeError(f"empty OUTCAR for spin task {task}: {outcar}")
        if text.count("Elapse") != 1 or text.count("TOTAL-FORCE") != 1:
            raise RuntimeError(
                f"spin task {task} did not finish normally; expected one "
                f"elapsed-time and one TOTAL-FORCE marker in: {outcar}"
            )
        oszicar = _require_file(stage / task / "OSZICAR", "spin OSZICAR", task)
        if not oszicar.read_text().strip():
            raise RuntimeError(f"empty OSZICAR for spin task {task}: {oszicar}")
    return len(paths)


def run_spin_tasks(jdata, mdata):
    """Submit saved tasks, retrieve OUTCAR/OSZICAR, and check termination."""
    make_spin_submission(jdata, mdata).run_submission()
    return check_spin_results(jdata)
