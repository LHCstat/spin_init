"""Generate VASP AIMD snapshots for the spin initialization workflow."""

import os
import shutil
import tempfile
from pathlib import Path

from dpgen import dlog
from dpgen.data.arginfo import spin_init_jdata_arginfo
from dpgen.data.gen import poscar_scale
from dpgen.data.tools.create_random_disturb import create_disturbs_ase_dev
from dpgen.dispatcher.Dispatcher import make_submission
from dpgen.generator.lib.utils import check_api_version
from dpgen.remote.decide_machine import convert_mdata
from dpgen.util import load_file, normalize

SCALE_PERT_DIR = "00.scale_pert"
MD_DIR = "01.md"
DISP_DIR = "02.disp"


def _expected_task_paths(jdata):
    """Return task paths in the stable scale/perturbation order."""
    return [
        f"scale-{scale:.3f}/{index:06d}"
        for scale in jdata["scale"]
        for index in range(jdata["pert_numb"] + 1)
    ]


def _unique(items):
    """Return values in first-seen order without duplicates."""
    return list(dict.fromkeys(items))


def validate_spin_init_parameters(jdata):
    """Validate workflow-specific ranges not expressed by the dargs schema."""
    if not jdata["stages"]:
        raise ValueError("stages must contain at least one stage")
    super_cell = jdata["super_cell"]
    if len(super_cell) != 3 or any(value <= 0 for value in super_cell):
        raise ValueError("super_cell must contain three positive integers")
    if not jdata["scale"] or any(value <= 0 for value in jdata["scale"]):
        raise ValueError("scale must contain at least one positive value")
    if jdata["pert_numb"] < 0:
        raise ValueError("pert_numb must be non-negative")
    if jdata["pert_box"] < 0:
        raise ValueError("pert_box must be non-negative")
    if jdata["pert_atom"] < 0:
        raise ValueError("pert_atom must be non-negative")
    if jdata["md_nstep"] < 0:
        raise ValueError("md_nstep must be non-negative")
    if not jdata["potcars"]:
        raise ValueError("potcars must contain at least one file")


def _require_file(path, description, task=None):
    """Return a resolved input file path or raise a contextual error."""
    file_path = Path(path).resolve()
    context = f" for task {task}" if task is not None else ""
    if not file_path.is_file():
        raise FileNotFoundError(f"{description} does not exist{context}: {file_path}")
    return file_path


def _make_relative_symlink(link_path, target_path, final_link_path=None):
    """Create a real relative symbolic link without a copy fallback."""
    if link_path.exists() or link_path.is_symlink():
        raise RuntimeError(f"symbolic link path already exists: {link_path}")
    link_reference = final_link_path if final_link_path is not None else link_path
    relative_target = os.path.relpath(target_path, start=link_reference.parent)
    os.symlink(relative_target, link_path)


def _synchronize_md_nstep(jdata):
    """Follow VASP INCAR NSW, matching the existing init_bulk behavior."""
    from pymatgen.io.vasp.inputs import Incar

    incar_path = _require_file(jdata["md_incar"], "MD INCAR source")
    incar = Incar.from_file(incar_path)
    if "NSW" in incar and int(incar["NSW"]) != jdata["md_nstep"]:
        dlog.warning(
            "spin_init md_nstep=%s differs from NSW=%s in %s; using NSW",
            jdata["md_nstep"],
            incar["NSW"],
            incar_path,
        )
        jdata["md_nstep"] = int(incar["NSW"])


def make_spin_init_structures(jdata):
    """Create scaled and perturbed POSCAR tasks from one input structure."""
    from pymatgen.core import Structure

    poscar_source = _require_file(jdata["from_poscar_path"], "POSCAR source")
    output_root = Path(jdata.get("out_dir", ".")).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    stage_path = output_root / SCALE_PERT_DIR
    if stage_path.exists():
        raise RuntimeError(f"output stage already exists: {stage_path}")
    task_paths = []

    with tempfile.TemporaryDirectory(dir=str(output_root)) as temporary_directory:
        scratch_stage = Path(temporary_directory) / SCALE_PERT_DIR
        scratch_stage.mkdir()
        supercell_path = Path(temporary_directory) / "POSCAR.supercell"
        structure = Structure.from_file(poscar_source)
        structure.make_supercell(jdata["super_cell"])
        structure.to(filename=str(supercell_path), fmt="poscar")

        for scale in jdata["scale"]:
            scale_name = f"scale-{scale:.3f}"
            scale_path = scratch_stage / scale_name
            scale_path.mkdir()
            scaled_poscar = scale_path / "POSCAR"
            poscar_scale(str(supercell_path), str(scaled_poscar), scale)

            previous_directory = Path.cwd()
            try:
                os.chdir(scale_path)
                create_disturbs_ase_dev(
                    "POSCAR",
                    jdata["pert_numb"],
                    dmax=jdata["pert_atom"],
                    etmax=jdata["pert_box"],
                    ofmt="vasp",
                )
            finally:
                os.chdir(previous_directory)

            unperturbed_path = scale_path / "000000"
            unperturbed_path.mkdir()
            shutil.move(str(scaled_poscar), str(unperturbed_path / "POSCAR"))
            task_paths.append(f"{scale_name}/000000")

            for index in range(1, jdata["pert_numb"] + 1):
                perturbation_path = scale_path / f"{index:06d}"
                perturbation_path.mkdir()
                shutil.move(
                    str(scale_path / f"POSCAR{index}.vasp"),
                    str(perturbation_path / "POSCAR"),
                )
                task_paths.append(f"{scale_name}/{index:06d}")

        scratch_stage.replace(stage_path)

    return task_paths


def make_spin_init_md(jdata, mdata):
    """Create VASP AIMD tasks with relative POSCAR/POTCAR symbolic links."""
    output_root = Path(jdata.get("out_dir", ".")).resolve()
    structure_stage = output_root / SCALE_PERT_DIR
    task_paths = _expected_task_paths(jdata)

    source_poscars = {}
    for task in task_paths:
        source_poscars[task] = _require_file(
            structure_stage / task / "POSCAR", "POSCAR source", task
        )
    md_incar = _require_file(jdata["md_incar"], "MD INCAR source")
    potcars = [
        _require_file(path, "POTCAR source") for path in jdata.get("potcars", [])
    ]
    if not potcars:
        raise ValueError("at least one POTCAR source is required")
    for path in mdata.get("fp_user_forward_files", []):
        _require_file(path, "user forward file")

    stage_path = output_root / MD_DIR
    if stage_path.exists():
        raise RuntimeError(f"output stage already exists: {stage_path}")

    output_root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=str(output_root)) as temporary_directory:
        scratch_stage = Path(temporary_directory) / MD_DIR
        scratch_stage.mkdir()
        common_incar = scratch_stage / "INCAR"
        shutil.copy2(md_incar, common_incar)
        common_potcar = scratch_stage / "POTCAR"
        with common_potcar.open("wb") as output_file:
            for potcar in potcars:
                with potcar.open("rb") as input_file:
                    shutil.copyfileobj(input_file, output_file)

        for task in task_paths:
            task_path = scratch_stage / task
            task_path.mkdir(parents=True)

            # These targets describe the final directory tree. Moving the complete
            # scratch stage into place preserves every relative link.
            final_task_path = stage_path / task
            _make_relative_symlink(
                task_path / "POSCAR",
                structure_stage / task / "POSCAR",
                final_task_path / "POSCAR",
            )
            _make_relative_symlink(
                task_path / "INCAR",
                stage_path / "INCAR",
                final_task_path / "INCAR",
            )
            _make_relative_symlink(
                task_path / "POTCAR",
                stage_path / "POTCAR",
                final_task_path / "POTCAR",
            )

        scratch_stage.replace(stage_path)

    # User-forwarded files (for example vasp.slurm) are task inputs only. The
    # command in machine.json independently decides whether they are executed.
    from dpgen.generator.lib.utils import symlink_user_forward_files

    symlink_user_forward_files(
        mdata=mdata,
        task_type="fp",
        work_path=str(stage_path),
        task_format={"fp": "scale-*/00*"},
    )
    return task_paths


def run_spin_init_md(jdata, mdata):
    """Submit all spin-init AIMD tasks through dpdispatcher."""
    work_path = Path(jdata.get("out_dir", ".")).resolve() / MD_DIR
    if not work_path.is_dir():
        raise FileNotFoundError(f"MD stage does not exist: {work_path}")

    task_directories = sorted(
        path for path in work_path.glob("scale-*/00*") if path.is_dir()
    )
    if not task_directories:
        raise RuntimeError(f"no AIMD tasks found under: {work_path}")
    run_tasks = [path.relative_to(work_path).as_posix() for path in task_directories]

    user_forward_files = mdata.get("fp_user_forward_files", [])
    forward_files = _unique(
        ["POSCAR", "INCAR", "POTCAR"]
        + [os.path.basename(path) for path in user_forward_files]
    )
    backward_files = _unique(
        ["OUTCAR", "XDATCAR"] + mdata.get("fp_user_backward_files", [])
    )

    check_api_version(mdata)
    submission = make_submission(
        mdata["fp_machine"],
        mdata["fp_resources"],
        commands=[mdata["fp_command"]],
        work_path=str(work_path),
        run_tasks=run_tasks,
        group_size=mdata["fp_group_size"],
        forward_common_files=[],
        forward_files=forward_files,
        backward_files=backward_files,
        outlog="fp.log",
        errlog="fp.log",
    )
    submission.run_submission()
    return run_tasks


def _task_context(task_path):
    """Return the scale/task portion of a task path for diagnostics."""
    path = Path(task_path)
    if len(path.parts) >= 2:
        return f"{path.parent.name}/{path.name}"
    return str(path)


def check_vasp_md_complete(task_path, md_nstep):
    """Check OUTCAR completion without inspecting or parsing XDATCAR."""
    task_path = Path(task_path).resolve()
    outcar = task_path / "OUTCAR"
    context = _task_context(task_path)
    if not outcar.is_file():
        raise FileNotFoundError(f"AIMD task {context} has no OUTCAR: {outcar}")
    text = outcar.read_text(errors="replace")
    if not text:
        raise RuntimeError(f"AIMD task {context} has an empty OUTCAR: {outcar}")
    if text.count("Elapse") != 1:
        raise RuntimeError(
            f"AIMD task {context} did not finish normally; expected one "
            f"elapsed-time marker in: {outcar}"
        )
    expected_force_blocks = 1 if md_nstep == 0 else md_nstep
    force_blocks = text.count("TOTAL-FORCE")
    if force_blocks != expected_force_blocks:
        raise RuntimeError(
            f"AIMD task {context} is incomplete: {outcar} contains "
            f"{force_blocks} TOTAL-FORCE blocks, expected {expected_force_blocks}"
        )


def parse_xdatcar_snapshots(xdatcar_path):
    """Parse all XDATCAR frames with pymatgen and validate frame consistency."""
    from pymatgen.io.vasp.outputs import Xdatcar

    xdatcar_path = Path(xdatcar_path).resolve()
    if not xdatcar_path.is_file():
        raise FileNotFoundError(f"XDATCAR does not exist: {xdatcar_path}")
    if xdatcar_path.stat().st_size == 0:
        raise RuntimeError(f"XDATCAR is empty: {xdatcar_path}")
    try:
        structures = list(Xdatcar(xdatcar_path).structures)
    except Exception as error:
        raise RuntimeError(
            f"failed to parse XDATCAR {xdatcar_path}: {error}"
        ) from error
    if not structures:
        raise RuntimeError(f"XDATCAR contains no trajectory frames: {xdatcar_path}")

    reference_species = [site.species_string for site in structures[0]]
    reference_count = len(reference_species)
    for frame_index, structure in enumerate(structures, start=1):
        species = [site.species_string for site in structure]
        if len(species) != reference_count or species != reference_species:
            raise RuntimeError(
                f"inconsistent frame {frame_index} in XDATCAR {xdatcar_path}: "
                f"expected {reference_count} atoms with species {reference_species}, "
                f"found {len(species)} atoms with species {species}"
            )
    return structures


def collect_xdatcar_snapshots(jdata):
    """Write every complete AIMD trajectory frame to ``02.disp``."""
    from pymatgen.io.vasp.inputs import Poscar

    output_root = Path(jdata.get("out_dir", ".")).resolve()
    md_stage = output_root / MD_DIR
    disp_stage = output_root / DISP_DIR
    if disp_stage.exists():
        raise RuntimeError(f"output stage already exists: {disp_stage}")
    if not md_stage.is_dir():
        raise FileNotFoundError(f"MD stage does not exist: {md_stage}")

    task_paths = _expected_task_paths(jdata)
    if not task_paths:
        raise RuntimeError(f"no AIMD tasks configured under: {md_stage}")

    # Complete every check before creating 02.disp so one bad task cannot leave
    # a partially collected dataset.
    trajectories = {}
    for task in task_paths:
        task_path = md_stage / task
        check_vasp_md_complete(task_path, jdata["md_nstep"])
        trajectories[task] = parse_xdatcar_snapshots(task_path / "XDATCAR")

    frame_count = 0
    with tempfile.TemporaryDirectory(dir=str(output_root)) as temporary_directory:
        scratch_stage = Path(temporary_directory) / DISP_DIR
        scratch_stage.mkdir()
        for task, structures in trajectories.items():
            for frame_index, structure in enumerate(structures):
                frame_path = scratch_stage / task / f"{frame_index:02d}"
                frame_path.mkdir(parents=True)
                Poscar(structure).write_file(frame_path / "POSCAR")
                frame_count += 1
        scratch_stage.replace(disp_stage)
    return frame_count


def gen_spin_init(args):
    """Run the spin initialization workflow."""
    jdata = normalize(spin_init_jdata_arginfo(), load_file(args.PARAM))
    validate_spin_init_parameters(jdata)
    stages = [int(stage) for stage in jdata["stages"]]
    for stage in stages:
        if stage not in (1, 2, 3):
            raise RuntimeError(f"unknown spin_init stage {stage}")
    _synchronize_md_nstep(jdata)

    mdata = None
    if args.MACHINE is not None:
        # Match init_bulk: convert either the current nested ``fp`` schema or
        # the legacy flat ``fp_*`` schema into the fields used by the runner.
        mdata = load_file(args.MACHINE)
        mdata = convert_mdata(mdata, ["fp"])

    output_root = Path(jdata["out_dir"]).resolve()
    output_root.mkdir(parents=True, exist_ok=True)
    parameter_source = Path(args.PARAM).resolve()
    parameter_copy = output_root / "param.json"
    if parameter_source != parameter_copy.resolve():
        shutil.copy2(parameter_source, parameter_copy)

    for stage in stages:
        if stage == 1:
            make_spin_init_structures(jdata)
        elif stage == 2:
            make_spin_init_md(jdata, mdata or {})
            if mdata is not None:
                run_spin_init_md(jdata, mdata)
        elif stage == 3:
            collect_xdatcar_snapshots(jdata)
