"""Stage 5: RMSE selection and scale-wise convert-data orchestration."""

import json
import os
import shutil
import tempfile
from pathlib import Path

from dpgen import dlog
from dpgen.data.out2npy import convert_extxyz_to_raw
from dpgen.data.spin_data import RMSE_LIMIT, magnetic_rmse, read_final_magnetization
from dpgen.data.spin_tasks import _saved_tasks, check_spin_results, read_spin_incar
from dpgen.dispatcher.Dispatcher import make_submission
from dpgen.generator.lib.utils import check_api_version
from dpgen.remote.decide_machine import convert_mdata


def select_spin_tasks(jdata):
    """Check completion, then retain tasks whose final moment RMSE passes."""
    from pymatgen.io.vasp.inputs import Poscar

    stage, tasks = _saved_tasks(jdata)
    check_spin_results(jdata)
    selected = []
    rejected = []
    next_index = {}
    for task in tasks:
        folder = stage / task
        try:
            natoms = len(Poscar.from_file(folder / "POSCAR").structure)
            _, initial = read_spin_incar(folder / "INCAR", natoms)
            _, final = read_final_magnetization(folder / "OUTCAR", natoms)
            rmse = magnetic_rmse(initial, final)
        except Exception as error:
            raise RuntimeError(f"spin data task {task}: {error}") from error
        if rmse > RMSE_LIMIT:
            rejected.append(
                {"task": task, "rmse": rmse, "reason": "RMSE exceeds 5.0e-3"}
            )
            dlog.warning("spin data task %s rejected: RMSE %.8g", task, rmse)
            continue
        scale = task.split("/", 1)[0]
        next_index[scale] = next_index.get(scale, 0) + 1
        selected.append(
            {
                "task": task,
                "scale": scale,
                "index": next_index[scale],
                "rmse": rmse,
            }
        )
    if not selected:
        raise RuntimeError(
            f"no spin tasks passed RMSE <= {RMSE_LIMIT:g}; "
            f"rejected {len(rejected)} of {len(tasks)} tasks"
        )
    return {"rmse_limit": RMSE_LIMIT, "selected": selected, "rejected": rejected}


def prepare_converter_inputs(stage, scratch, selection):
    """Make upload-safe data/OUTCAR-N and data/OSZICAR-N input pairs."""
    stage = Path(stage)
    scratch = Path(scratch)
    scales = []
    for row in selection["selected"]:
        scale = row["scale"]
        if scale not in scales:
            scales.append(scale)
        source = stage / row["task"]
        data = scratch / scale / "data"
        data.mkdir(parents=True, exist_ok=True)
        for name in ("OUTCAR", "OSZICAR"):
            original = source / name
            if not original.is_file() or original.stat().st_size == 0:
                raise FileNotFoundError(
                    f"spin data task {row['task']}: missing or empty {original}"
                )
            shutil.copy2(original, data / f"{name}-{row['index']}")
    (scratch / "selection.json").write_text(
        json.dumps(selection, indent=2), encoding="utf-8"
    )
    return scales


def run_convert_data(work_path, scales, mdata):
    """Submit one ``nequip-data`` task per scale through dpdispatcher."""
    if not mdata or "convert-data" not in mdata:
        raise ValueError("Stage 5 requires MACHINE with a convert-data entry")
    entry = mdata["convert-data"]
    if isinstance(entry, (list, tuple)):
        if not entry:
            raise ValueError("convert-data requires a nonempty machine entry")
        entry = entry[0]
    if not isinstance(entry, dict):
        raise ValueError("convert-data must be a machine object or nonempty list")
    if not isinstance(entry.get("command"), str) or not entry["command"].strip():
        raise ValueError("convert-data requires a nonempty command")
    for field in ("machine", "resources"):
        if not isinstance(entry.get(field), dict):
            raise ValueError(f"convert-data requires a {field} object")
    converted = convert_mdata(dict(mdata), ["convert-data"])
    prefix = "convert-data_"
    check_api_version(converted)
    context = converted[prefix + "machine"].get("context_type", "").lower()
    if context in (
        "lebesguecontext",
        "dpcloudservercontext",
        "openapicontext",
        "bohriumcontext",
    ):
        try:
            import oss2  # noqa: F401
        except ImportError as error:
            raise RuntimeError(
                'Bohrium dependency unavailable; run: python -m pip install "dpdispatcher[bohrium]"'
            ) from error
    # Some dpdispatcher contexts move the nested backward file into the local
    # task directory without creating its parent first.
    for scale in scales:
        (Path(work_path) / scale / "out").mkdir(parents=True, exist_ok=True)
    submission = make_submission(
        converted[prefix + "machine"],
        converted[prefix + "resources"],
        commands=[converted[prefix + "command"]],
        work_path=str(work_path),
        run_tasks=scales,
        group_size=converted[prefix + "group_size"],
        forward_common_files=[],
        forward_files=["data"],
        backward_files=["out/data.extxyz"],
        outlog="convert.log",
        errlog="convert.log",
    )
    submission.run_submission()


def convert_scale_outputs(work_path, scales):
    """Place extxyz, raw, and npy together in each PDF-specified out/."""
    work_path = Path(work_path)
    counts = {}
    for scale in scales:
        out = work_path / scale / "out"
        extxyz = out / "data.extxyz"
        if not extxyz.is_file() or extxyz.stat().st_size == 0:
            raise FileNotFoundError(
                f"convert-data did not return a nonempty file for {scale}: {extxyz}"
            )
        counts[scale] = convert_extxyz_to_raw(extxyz, out)
    return counts


def publish_scale_links(stage, scratch, destination, selection):
    """Replace upload copies with the PDF's two levels of relative links."""
    stage = Path(stage)
    scratch = Path(scratch)
    destination = Path(destination)
    scales = list(dict.fromkeys(row["scale"] for row in selection["selected"]))
    for scale in scales:
        final_data = stage / scale / "data"
        if final_data.exists() or final_data.is_symlink():
            raise FileExistsError(f"spin data directory already exists: {final_data}")
    created = []
    try:
        for scale in scales:
            final_data = stage / scale / "data"
            final_data.mkdir()
            created.append(final_data)
            for row in selection["selected"]:
                if row["scale"] != scale:
                    continue
                for name in ("OUTCAR", "OSZICAR"):
                    source = stage / row["task"] / name
                    target = final_data / f"{name}-{row['index']}"
                    os.symlink(os.path.relpath(source, start=final_data), target)
            upload_data = scratch / scale / "data"
            shutil.rmtree(upload_data)
            final_link = destination / scale / "data"
            os.symlink(
                os.path.relpath(final_data, start=final_link.parent),
                upload_data,
                target_is_directory=True,
            )
    except Exception:
        for path in reversed(created):
            shutil.rmtree(path)
        raise


def collect_spin_data(jdata, mdata):
    """Run RMSE selection, convert-data, and extxyz/raw/npy export."""
    root = Path(jdata.get("out_dir", ".")).resolve()
    destination = root / "04.data"
    if destination.exists() or destination.is_symlink():
        raise RuntimeError(f"output stage already exists: {destination}")
    if not mdata or "convert-data" not in mdata:
        raise ValueError("Stage 5 requires MACHINE with a convert-data entry")
    stage, _ = _saved_tasks(jdata)
    selection = select_spin_tasks(jdata)
    scales = list(dict.fromkeys(row["scale"] for row in selection["selected"]))
    for scale in scales:
        data = stage / scale / "data"
        if data.exists() or data.is_symlink():
            raise FileExistsError(f"spin data directory already exists: {data}")
    with tempfile.TemporaryDirectory(dir=root) as temporary:
        scratch = Path(temporary) / "04.data"
        scratch.mkdir()
        prepare_converter_inputs(stage, scratch, selection)
        run_convert_data(scratch, scales, mdata)
        counts = convert_scale_outputs(scratch, scales)
        publish_scale_links(stage, scratch, destination, selection)
        try:
            scratch.replace(destination)
        except Exception:
            for scale in scales:
                shutil.rmtree(stage / scale / "data")
            raise
    return sum(counts.values())
