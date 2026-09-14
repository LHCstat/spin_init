"""Magnetic-moment perturbations for ``spin_init`` Stage 4."""

from itertools import product

import numpy as np

_CARTESIAN_AXES = np.eye(3)
_COLLINEAR_TOLERANCE = 8 * np.finfo(float).eps


def _unit_and_magnitude(vector):
    """Return a scale-safe unit vector and magnitude for a finite 3-vector."""
    scale = float(np.max(np.abs(vector)))
    if scale == 0.0:
        return None, 0.0
    scaled = vector / scale
    scaled_norm = float(np.sqrt(np.dot(scaled, scaled)))
    unit = scaled / scaled_norm
    if scale > np.finfo(float).max / scaled_norm:
        magnitude = float("inf")
    else:
        magnitude = scale * scaled_norm
    return unit, magnitude


def _moments_array(moments):
    try:
        result = np.asarray(moments, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError("magnetic moments must be numeric") from error
    if result.ndim != 2 or result.shape[1] != 3 or not np.isfinite(result).all():
        raise ValueError("magnetic moments must be a finite N x 3 array")
    return result.copy()


def _fallback_direction(unit_moment):
    index = int(np.argmin(np.abs(_CARTESIAN_AXES @ unit_moment)))
    return _CARTESIAN_AXES[index]


def cant_moments(moments, angle, rcut, direction=None):
    """Cant every nonzero moment while preserving its original magnitude."""
    result = _moments_array(moments)
    if isinstance(angle, (bool, np.bool_)) or isinstance(rcut, (bool, np.bool_)):
        raise ValueError("Canting angle and Rcut must be numeric, not boolean")
    try:
        angle = float(angle)
        rcut = float(rcut)
    except (TypeError, ValueError) as error:
        raise ValueError("Canting angle and Rcut must be numeric") from error
    if not np.isfinite(angle):
        raise ValueError("Canting angle must be finite")
    if not np.isfinite(rcut) or rcut < 0:
        raise ValueError("Canting Rcut must be a finite non-negative number")

    supplied_direction = None
    if direction is not None:
        try:
            supplied_direction = np.asarray(direction, dtype=float)
        except (TypeError, ValueError) as error:
            raise ValueError("Canting direction must be numeric") from error
        if (
            supplied_direction.shape != (3,)
            or not np.isfinite(supplied_direction).all()
        ):
            raise ValueError("Canting direction must be one finite 3-vector")
        supplied_direction, _ = _unit_and_magnitude(supplied_direction)

    radians = np.deg2rad(angle)
    cosine = np.cos(radians)
    sine = np.sin(radians)
    for atom_index, moment in enumerate(result):
        unit_moment, magnitude = _unit_and_magnitude(moment)
        if unit_moment is None:
            continue
        if not np.isfinite(magnitude):
            raise ValueError(
                f"atom {atom_index}: magnetic-moment magnitude is too large"
            )
        if rcut > magnitude:
            raise ValueError(
                f"atom {atom_index}: Rcut={rcut:g} exceeds magnetic-moment "
                f"magnitude {magnitude:g}"
            )

        reference = supplied_direction
        if reference is None:
            reference = _fallback_direction(unit_moment)
        projection = reference - np.dot(reference, unit_moment) * unit_moment
        unit_projection, projection_magnitude = _unit_and_magnitude(projection)
        if projection_magnitude <= _COLLINEAR_TOLERANCE:
            reference = _fallback_direction(unit_moment)
            projection = reference - np.dot(reference, unit_moment) * unit_moment
            unit_projection, _ = _unit_and_magnitude(projection)

        transverse_direction = cosine * unit_projection + sine * np.cross(
            unit_moment, unit_projection
        )
        transverse_magnitude = rcut
        ratio = transverse_magnitude / magnitude
        longitudinal = magnitude * np.sqrt(max(0.0, (1.0 - ratio) * (1.0 + ratio)))
        result[atom_index] = (
            longitudinal * unit_moment + transverse_magnitude * transverse_direction
        )
    return result


def _number_values(value, field, non_negative=False):
    objects = np.asarray(value, dtype=object)
    if any(isinstance(item, (bool, np.bool_)) for item in objects.flat):
        raise ValueError(f"Canting {field} values must not be boolean")
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"Canting {field} must be numeric") from error
    if array.ndim == 0:
        array = array.reshape(1)
    if array.ndim != 1 or not array.size or not np.isfinite(array).all():
        raise ValueError(f"Canting {field} must be a finite number or nonempty list")
    if non_negative and np.any(array < 0):
        raise ValueError(f"Canting {field} values must be non-negative")
    return [float(item) for item in array]


def _direction_values(value):
    if value is None:
        return [None]
    objects = np.asarray(value, dtype=object)
    if any(isinstance(item, (bool, np.bool_)) for item in objects.flat):
        raise ValueError("Canting direction values must not be boolean")
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError("Canting direction must be numeric") from error
    if array.shape == (3,):
        array = array.reshape(1, 3)
    if array.ndim != 2 or array.shape[1] != 3 or not len(array):
        raise ValueError("Canting direction must be a 3-vector or a list of 3-vectors")
    if not np.isfinite(array).all():
        raise ValueError("Canting direction values must be finite")
    return [row.copy() for row in array]


def build_spin_perturbation(pert_spin):
    """Validate Canting blocks and return ``(count, provider)`` for Stage 4."""
    if pert_spin is None:
        return 0, None
    if not isinstance(pert_spin, list):
        raise ValueError("pert_spin must be a list of perturbation blocks")
    if not pert_spin:
        return 0, None

    combinations = []
    for block_index, block in enumerate(pert_spin):
        if not isinstance(block, dict) or not block:
            raise ValueError(f"pert_spin[{block_index}] must be a nonempty object")
        for mode, parameters in block.items():
            if mode != "Canting":
                raise NotImplementedError(
                    f"spin perturbation mode {mode!r} is not implemented"
                )
            if not isinstance(parameters, dict):
                raise ValueError(f"pert_spin[{block_index}].Canting must be an object")
            unknown = set(parameters) - {"angle", "Rcut", "direction"}
            if unknown:
                raise ValueError(
                    f"unknown Canting parameter(s): {', '.join(sorted(unknown))}"
                )
            missing = {"angle", "Rcut"} - set(parameters)
            if missing:
                raise ValueError(
                    f"missing Canting parameter(s): {', '.join(sorted(missing))}"
                )
            angles = _number_values(parameters["angle"], "angle")
            rcuts = _number_values(parameters["Rcut"], "Rcut", non_negative=True)
            directions = _direction_values(parameters.get("direction"))
            combinations.extend(product(angles, rcuts, directions))

    def provider(moments, count):
        if count != len(combinations):
            raise ValueError(
                f"Canting provider expected {len(combinations)} configurations, "
                f"received {count}"
            )
        configurations = {}
        for index, (angle, rcut, direction) in enumerate(combinations, start=1):
            name = f"C{index}"
            try:
                configurations[name] = cant_moments(
                    moments, angle=angle, rcut=rcut, direction=direction
                )
            except ValueError as error:
                raise ValueError(f"{name}: {error}") from error
        return configurations

    return len(combinations), provider
