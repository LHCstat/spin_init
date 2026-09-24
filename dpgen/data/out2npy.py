"""Convert multi-frame magnetic extxyz to raw files and ``set.000/*.npy``."""

import argparse
import shlex
import shutil
import tempfile
from contextlib import ExitStack
from pathlib import Path

import numpy as np


FRAME_FILES = (
    "box.raw",
    "coord.raw",
    "energy.raw",
    "force.raw",
    "force_mag.raw",
    "spin.raw",
    "virial.raw",
)
REQUIRED_PROPERTIES = {
    "species": 1,
    "pos": 3,
    "spin_length": 1,
    "initial_magmoms": 3,
    "spin_forces_vert": 3,
    "forces": 3,
}


def _numbers(values, context):
    try:
        result = np.asarray([float(value) for value in values], dtype=float)
    except ValueError as error:
        raise ValueError(f"{context}: invalid numeric value") from error
    if not np.isfinite(result).all():
        raise ValueError(f"{context}: non-finite numeric value")
    return result


def _header(line, context):
    try:
        fields = dict(token.split("=", 1) for token in shlex.split(line))
        properties = fields["Properties"].split(":")
        if len(properties) % 3:
            raise ValueError("invalid Properties triplets")
        columns = {}
        offset = 0
        for index in range(0, len(properties), 3):
            name, kind, width_text = properties[index : index + 3]
            width = int(width_text)
            if name in columns or width < 1 or kind not in ("S", "R", "I"):
                raise ValueError("invalid Properties definition")
            columns[name] = (offset, width)
            offset += width
        for name, width in REQUIRED_PROPERTIES.items():
            if name not in columns or columns[name][1] != width:
                raise ValueError(f"Properties requires {name} with width {width}")
        cell = _numbers(fields["Lattice"].split(), context + " Lattice")
        stress = _numbers(fields["stress"].split(), context + " stress")
        energy = _numbers([fields["energy"]], context + " energy")
        if cell.size != 9 or stress.size != 9:
            raise ValueError("Lattice and stress each require nine values")
        cell = cell.reshape(3, 3)
        stress = stress.reshape(3, 3)
        # The input is always a 3x3 VASP cell.  Evaluate its triple product
        # directly, without requiring NumPy's native linear-algebra backend.
        a, b, c = cell
        volume = abs(
            float(
                a[0] * (b[1] * c[2] - b[2] * c[1])
                - a[1] * (b[0] * c[2] - b[2] * c[0])
                + a[2] * (b[0] * c[1] - b[1] * c[0])
            )
        )
        if not np.isfinite(volume) or volume <= 0:
            raise ValueError("Lattice has zero or invalid volume")
    except (KeyError, TypeError, ValueError) as error:
        raise ValueError(f"{context}: invalid extxyz header: {error}") from error
    return columns, offset, cell, stress, energy, volume


def _read_frame(source, frame_number):
    count_line = source.readline()
    while count_line and not count_line.strip():
        count_line = source.readline()
    if not count_line:
        return None
    context = f"{source.name} frame {frame_number}"
    try:
        natoms = int(count_line.strip())
    except ValueError as error:
        raise ValueError(
            f"{context}: invalid atom count {count_line.strip()!r}"
        ) from error
    if natoms <= 0:
        raise ValueError(f"{context}: atom count must be positive")
    comment = source.readline()
    if not comment:
        raise ValueError(f"{context}: missing extxyz header")
    columns, width, cell, stress, energy, volume = _header(comment, context)
    species = []
    arrays = {name: [] for name in REQUIRED_PROPERTIES if name != "species"}
    for atom in range(natoms):
        row = source.readline()
        if not row:
            raise ValueError(f"{context}: incomplete frame at atom {atom + 1}/{natoms}")
        tokens = row.split()
        if len(tokens) != width:
            raise ValueError(
                f"{context}: atom {atom + 1} has {len(tokens)} columns, expected {width}"
            )
        start, size = columns["species"]
        species.append(tokens[start])
        for name in arrays:
            start, size = columns[name]
            arrays[name].append(
                _numbers(tokens[start : start + size], context + " " + name)
            )
    for name in arrays:
        arrays[name] = np.asarray(arrays[name], dtype=float)
    rows = {
        "box.raw": cell.reshape(-1),
        "coord.raw": arrays["pos"].reshape(-1),
        "energy.raw": energy,
        "force.raw": arrays["forces"].reshape(-1),
        "force_mag.raw": arrays["spin_forces_vert"].reshape(-1),
        "spin.raw": (arrays["spin_length"] * arrays["initial_magmoms"]).reshape(-1),
        "virial.raw": (-volume * stress).reshape(-1),
    }
    return species, rows


def _write_npy_files(scratch, frame_number):
    """Create float64 arrays from the raw rows without dropping single frames."""
    set_dir = scratch / "set.000"
    set_dir.mkdir()
    for name in FRAME_FILES:
        values = np.loadtxt(scratch / name, dtype=np.float64, ndmin=2)
        if values.shape[0] != frame_number:
            raise ValueError(
                f"{scratch / name}: {values.shape[0]} frames, expected {frame_number}"
            )
        array = values[:, 0] if name == "energy.raw" else values
        np.save(set_dir / name.replace(".raw", ".npy"), array)


def convert_extxyz_to_raw(input_path, output_dir):
    """Write raw and npy files, preserving frame and atom order."""
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    if not input_path.is_file():
        raise FileNotFoundError(f"extxyz input does not exist: {input_path}")
    if output_dir.is_symlink() or (output_dir.exists() and not output_dir.is_dir()):
        raise FileExistsError(f"output is not a regular directory: {output_dir}")
    artifacts = ["type_map.raw", "type.raw", *FRAME_FILES, "set.000"]
    if output_dir.exists():
        for name in artifacts:
            target = output_dir / name
            if target.exists() or target.is_symlink():
                raise FileExistsError(f"conversion output already exists: {target}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=output_dir.parent) as temporary:
        scratch = Path(temporary) / "raw"
        scratch.mkdir()
        with input_path.open(encoding="utf-8") as source, ExitStack() as stack:
            writers = {
                name: stack.enter_context((scratch / name).open("w", encoding="utf-8"))
                for name in FRAME_FILES
            }
            first_species = None
            frame_number = 0
            while True:
                frame = _read_frame(source, frame_number)
                if frame is None:
                    break
                species, rows = frame
                if first_species is None:
                    first_species = species
                elif species != first_species:
                    raise ValueError(
                        f"{input_path} frame {frame_number}: atom count/order/species "
                        "differs from the first frame"
                    )
                for name, values in rows.items():
                    np.savetxt(writers[name], values.reshape(1, -1), fmt="%.18e")
                frame_number += 1
        if not frame_number:
            raise ValueError(f"extxyz input is empty: {input_path}")
        atom_names = list(dict.fromkeys(first_species))
        atom_types = [atom_names.index(name) for name in first_species]
        (scratch / "type_map.raw").write_text(
            "\n".join(atom_names) + "\n", encoding="utf-8"
        )
        np.savetxt(scratch / "type.raw", atom_types, fmt="%d")
        _write_npy_files(scratch, frame_number)
        if output_dir.exists():
            published = []
            try:
                for name in artifacts:
                    source = scratch / name
                    target = output_dir / name
                    if target.exists() or target.is_symlink():
                        raise FileExistsError(f"conversion output already exists: {target}")
                    source.replace(target)
                    published.append(target)
            except Exception:
                for target in reversed(published):
                    if target.is_dir():
                        shutil.rmtree(target)
                    else:
                        target.unlink()
                raise
        else:
            scratch.replace(output_dir)
    return frame_number


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("extxyz", help="multi-frame data.extxyz")
    parser.add_argument("output", help="directory for raw files and set.000/*.npy")
    args = parser.parse_args(argv)
    try:
        convert_extxyz_to_raw(args.extxyz, args.output)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
