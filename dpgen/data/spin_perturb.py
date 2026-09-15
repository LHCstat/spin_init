"""Magnetic-moment perturbations for ``spin_init`` Stage 4."""

import numpy as np

_CARTESIAN_AXES = np.eye(3)


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


def cant_moments(moments, angle, rng=None):
    """Cant every nonzero moment by ``angle`` with random atomwise azimuths."""
    result = _moments_array(moments)
    if np.asarray(angle, dtype=object).ndim != 0:
        raise ValueError("cant_moments angle must be one scalar value")
    angle = _angle_values(angle)[0]
    if angle == 0.0:
        return result
    if angle == 180.0:
        nonzero = np.any(result != 0.0, axis=1)
        result[nonzero] = -result[nonzero]
        return result
    if rng is None:
        rng = np.random.default_rng()
    if not callable(getattr(rng, "uniform", None)):
        raise ValueError("Canting rng must provide uniform(low, high)")

    radians = np.deg2rad(angle)
    polar_cosine = np.cos(radians)
    polar_sine = np.sin(radians)
    for atom_index, moment in enumerate(result):
        unit_moment, magnitude = _unit_and_magnitude(moment)
        if unit_moment is None:
            continue
        if not np.isfinite(magnitude):
            raise ValueError(
                f"atom {atom_index}: magnetic-moment magnitude is too large"
            )
        reference = _fallback_direction(unit_moment)
        projection = reference - np.dot(reference, unit_moment) * unit_moment
        first_axis, _ = _unit_and_magnitude(projection)
        second_axis = np.cross(unit_moment, first_axis)
        azimuth = rng.uniform(0.0, 2.0 * np.pi)
        transverse_direction = (
            np.cos(azimuth) * first_axis + np.sin(azimuth) * second_axis
        )
        result[atom_index] = magnitude * (
            polar_cosine * unit_moment + polar_sine * transverse_direction
        )
    return result


def _angle_values(value):
    objects = np.asarray(value, dtype=object)
    if any(isinstance(item, (bool, np.bool_)) for item in objects.flat):
        raise ValueError("Canting angle values must not be boolean")
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError("Canting angle must be numeric") from error
    if array.ndim == 0:
        array = array.reshape(1)
    if array.ndim != 1 or not array.size or not np.isfinite(array).all():
        raise ValueError("Canting angle must be a finite number or nonempty list")
    if np.any(array < 0) or np.any(array > 180):
        raise ValueError("Canting angle values must be between 0 and 180 degrees")
    return [float(item) for item in array]


def _seed_value(value):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError("Canting seed must be a non-negative integer")
    if value < 0:
        raise ValueError("Canting seed must be a non-negative integer")
    return int(value)


def build_spin_perturbation(pert_spin):
    """Validate Canting blocks and return ``(count, provider)`` for Stage 4."""
    if pert_spin is None:
        return 0, None
    if not isinstance(pert_spin, list):
        raise ValueError("pert_spin must be a list of perturbation blocks")
    if not pert_spin:
        return 0, None

    configurations = []
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
            unknown = set(parameters) - {"angle", "seed"}
            if unknown:
                raise ValueError(
                    f"unknown Canting parameter(s): {', '.join(sorted(unknown))}"
                )
            missing = {"angle"} - set(parameters)
            if missing:
                raise ValueError(
                    f"missing Canting parameter(s): {', '.join(sorted(missing))}"
                )
            angles = _angle_values(parameters["angle"])
            seed = _seed_value(parameters["seed"]) if "seed" in parameters else None
            rng = np.random.default_rng(seed)
            configurations.extend((angle, rng) for angle in angles)

    def provider(moments, count):
        if count != len(configurations):
            raise ValueError(
                f"Canting provider expected {len(configurations)} configurations, "
                f"received {count}"
            )
        perturbed = {}
        for index, (angle, rng) in enumerate(configurations, start=1):
            name = f"C{index}"
            try:
                perturbed[name] = cant_moments(moments, angle=angle, rng=rng)
            except ValueError as error:
                raise ValueError(f"{name}: {error}") from error
        return perturbed

    return len(configurations), provider
