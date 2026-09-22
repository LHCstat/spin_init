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


def _number_values(value, label, minimum, maximum):
    objects = np.asarray(value, dtype=object)
    if any(isinstance(item, (bool, np.bool_)) for item in objects.flat):
        raise ValueError(f"{label} values must not be boolean")
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if array.ndim == 0:
        array = array.reshape(1)
    if array.ndim != 1 or not array.size or not np.isfinite(array).all():
        raise ValueError(f"{label} must be a finite number or nonempty list")
    if np.any(array < minimum) or np.any(array > maximum):
        raise ValueError(
            f"{label} values must be between {minimum:g} and {maximum:g} degrees"
        )
    return [float(item) for item in array]


def _axis_values(value, label):
    objects = np.asarray(value, dtype=object)
    if any(isinstance(item, (bool, np.bool_)) for item in objects.flat):
        raise ValueError(f"{label} values must not be boolean")
    try:
        array = np.asarray(value, dtype=float)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if array.ndim == 1:
        array = array.reshape(1, -1)
    if (
        array.ndim != 2
        or array.shape[0] == 0
        or array.shape[1] != 3
        or not np.isfinite(array).all()
    ):
        raise ValueError(f"{label} must be a finite 3-vector or nonempty list")
    units = []
    for index, axis in enumerate(array):
        unit, _ = _unit_and_magnitude(axis)
        if unit is None:
            raise ValueError(f"{label}[{index}] must be nonzero")
        units.append(unit)
    return units


def _single_finite_number(value, label):
    if isinstance(value, (bool, np.bool_)):
        raise ValueError(f"{label} must not be boolean")
    try:
        result = float(value)
    except (TypeError, ValueError) as error:
        raise ValueError(f"{label} must be numeric") from error
    if np.asarray(value, dtype=object).ndim != 0 or not np.isfinite(result):
        raise ValueError(f"{label} must be one finite number")
    return result


def _scale_deltas(pert, pert_step, location="Scale"):
    pert_value = _single_finite_number(pert, f"{location}.pert")
    step_value = _single_finite_number(pert_step, f"{location}.pert_step")
    if not 0.0 < pert_value < 1.0:
        raise ValueError(f"{location}.pert must satisfy 0 < pert < 1")
    if step_value <= 0.0:
        raise ValueError(f"{location}.pert_step must be positive")
    steps = int(round(pert_value / step_value))
    if steps < 1 or not np.isclose(
        steps * step_value, pert_value, rtol=1e-12, atol=1e-15
    ):
        raise ValueError(f"{location}.pert / pert_step must be an integer")
    negative = [-step_value * index for index in range(steps, 0, -1)]
    positive = [step_value * index for index in range(1, steps + 1)]
    return [float(value) for value in negative + positive]


def scale_moments(moments, delta):
    """Scale every moment by the relative factor ``1 + delta``."""
    delta_value = _single_finite_number(delta, "Scale delta")
    result = _moments_array(moments) * (1.0 + delta_value)
    if not np.isfinite(result).all():
        raise ValueError("scaled magnetic moments must remain finite")
    return result


def rotate_moments(moments, angle, axis):
    """Rotate moments around a global axis using the right-hand rule."""
    result = _moments_array(moments)
    if np.asarray(angle, dtype=object).ndim != 0:
        raise ValueError("rotate_moments angle must be one scalar value")
    degrees = _number_values(angle, "Rotation angle", 0.0, 360.0)[0]
    if np.asarray(axis, dtype=object).ndim != 1:
        raise ValueError("rotate_moments axis must be one 3-vector")
    axis_unit = _axis_values(axis, "Rotation axis")[0]
    if degrees == 0.0 or degrees == 360.0:
        return result
    radians = np.deg2rad(degrees)
    cosine = np.cos(radians)
    sine = np.sin(radians)
    return (
        result * cosine
        + np.cross(axis_unit, result) * sine
        + np.outer(result @ axis_unit, axis_unit) * (1.0 - cosine)
    )


def _require_uniform_rng(rng, mode):
    if not callable(getattr(rng, "uniform", None)):
        raise ValueError(f"{mode} rng must provide uniform(low, high)")


def randomize_moments(moments, rng=None):
    """Give every nonzero moment an independent uniform sphere direction."""
    result = _moments_array(moments)
    if rng is None:
        rng = np.random.default_rng()
    _require_uniform_rng(rng, "Random")
    for atom_index, moment in enumerate(result):
        _, magnitude = _unit_and_magnitude(moment)
        if magnitude == 0.0:
            continue
        if not np.isfinite(magnitude):
            raise ValueError(
                f"atom {atom_index}: magnetic-moment magnitude is too large"
            )
        z = rng.uniform(-1.0, 1.0)
        azimuth = rng.uniform(0.0, 2.0 * np.pi)
        radius = np.sqrt(max(0.0, 1.0 - z * z))
        result[atom_index] = magnitude * np.array(
            [radius * np.cos(azimuth), radius * np.sin(azimuth), z]
        )
    return result


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
    _require_uniform_rng(rng, "Canting")

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


def rota_cant_moments(moments, rotation_angle, axis, canting_angle, rng=None):
    """Apply a global Rotation followed by atomwise Canting."""
    rotated = rotate_moments(moments, rotation_angle, axis)
    return cant_moments(rotated, canting_angle, rng=rng)


def _angle_values(value):
    return _number_values(value, "Canting angle", 0.0, 180.0)


def _seed_value(value, label="Canting seed"):
    if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
        raise ValueError(f"{label} must be a non-negative integer")
    if value < 0:
        raise ValueError(f"{label} must be a non-negative integer")
    return int(value)


def _validate_parameters(parameters, allowed, required, location):
    if not isinstance(parameters, dict):
        raise ValueError(f"{location} must be an object")
    unknown = set(parameters) - set(allowed)
    if unknown:
        raise ValueError(
            f"{location}: unknown parameter(s): {', '.join(sorted(unknown))}"
        )
    missing = set(required) - set(parameters)
    if missing:
        raise ValueError(
            f"{location}: missing parameter(s): {', '.join(sorted(missing))}"
        )


def _canting_variants(parameters, location):
    _validate_parameters(parameters, {"angle", "seed"}, {"angle"}, location)
    angles = _number_values(parameters["angle"], f"{location}.angle", 0.0, 180.0)
    seed = (
        _seed_value(parameters["seed"], f"{location}.seed")
        if "seed" in parameters
        else None
    )
    rng = np.random.default_rng(seed)
    return [
        (
            f"C{index}",
            lambda moments, angle=angle, rng=rng: cant_moments(
                moments, angle=angle, rng=rng
            ),
        )
        for index, angle in enumerate(angles, start=1)
    ]


def _rotation_variants(parameters, location):
    _validate_parameters(parameters, {"angle", "axis"}, {"angle", "axis"}, location)
    angles = _number_values(parameters["angle"], f"{location}.angle", 0.0, 360.0)
    axes = _axis_values(parameters["axis"], f"{location}.axis")
    variants = []
    for angle in angles:
        for axis in axes:
            name = f"R{len(variants) + 1}"
            variants.append(
                (
                    name,
                    lambda moments, angle=angle, axis=axis: rotate_moments(
                        moments, angle=angle, axis=axis
                    ),
                )
            )
    return variants


def _random_variants(parameters, location):
    _validate_parameters(parameters, {"num", "seed"}, {"num"}, location)
    number = parameters["num"]
    if (
        isinstance(number, (bool, np.bool_))
        or not isinstance(number, (int, np.integer))
        or number <= 0
    ):
        raise ValueError(f"{location}.num must be a positive integer")
    seed = (
        _seed_value(parameters["seed"], f"{location}.seed")
        if "seed" in parameters
        else None
    )
    rng = np.random.default_rng(seed)
    return [
        (
            f"Rand{index}",
            lambda moments, rng=rng: randomize_moments(moments, rng),
        )
        for index in range(1, int(number) + 1)
    ]


def _scale_variants(parameters, location):
    _validate_parameters(
        parameters, {"pert", "pert_step"}, {"pert", "pert_step"}, location
    )
    deltas = _scale_deltas(
        parameters["pert"], parameters["pert_step"], location=location
    )
    return [
        (
            f"S{index}",
            lambda moments, delta=delta: scale_moments(moments, delta),
        )
        for index, delta in enumerate(deltas, start=1)
    ]


def _rota_cant_variants(parameters, location):
    _validate_parameters(
        parameters,
        {"R_angle", "axis", "C_angle", "seed"},
        {"R_angle", "axis", "C_angle"},
        location,
    )
    rotation_angles = _number_values(
        parameters["R_angle"], f"{location}.R_angle", 0.0, 360.0
    )
    axes = _axis_values(parameters["axis"], f"{location}.axis")
    canting_angles = _number_values(
        parameters["C_angle"], f"{location}.C_angle", 0.0, 180.0
    )
    seed = (
        _seed_value(parameters["seed"], f"{location}.seed")
        if "seed" in parameters
        else None
    )
    rng = np.random.default_rng(seed)
    variants = []
    for rotation_angle in rotation_angles:
        for axis in axes:
            for canting_angle in canting_angles:
                name = f"RC{len(variants) + 1}"
                variants.append(
                    (
                        name,
                        lambda moments, rotation_angle=rotation_angle, axis=axis, canting_angle=canting_angle, rng=rng: (
                            rota_cant_moments(
                                moments,
                                rotation_angle,
                                axis,
                                canting_angle,
                                rng,
                            )
                        ),
                    )
                )
    return variants


_OPERATION_COMPILERS = {
    "Rotation": _rotation_variants,
    "Canting": _canting_variants,
    "Rota_Cant": _rota_cant_variants,
    "Random": _random_variants,
    "Scale": _scale_variants,
}


def _compile_operations(pert_spin):
    if pert_spin is None:
        return []
    if not isinstance(pert_spin, list):
        raise ValueError("pert_spin must be a list of perturbation blocks")
    operations = []
    for block_index, block in enumerate(pert_spin):
        block_location = f"pert_spin[{block_index}]"
        if not isinstance(block, dict) or not block:
            raise ValueError(f"{block_location} must be a nonempty object")
        if len(block) != 1:
            raise ValueError(f"{block_location} must contain exactly one operation")
        mode, parameters = next(iter(block.items()))
        location = f"{block_location}.{mode}"
        try:
            compiler = _OPERATION_COMPILERS[mode]
        except KeyError as error:
            raise NotImplementedError(
                f"{block_location}: spin perturbation mode {mode!r} is not implemented"
            ) from error
        group = f"{block_index:03d}-{mode.lower()}"
        operations.append((group, compiler(parameters, location)))
    return operations


def build_spin_perturbation(pert_spin):
    """Compile independent operation groups into the Stage 4 provider contract.

    Each block starts from the input moments. Only parameter variants inside
    that block form a Cartesian product; different blocks never compose.
    """
    operations = _compile_operations(pert_spin)
    if not operations:
        return 0, None
    count = sum(len(variants) for _, variants in operations)

    def provider(moments, requested_count):
        if requested_count != count:
            raise ValueError(
                f"spin provider expected {count} configurations, "
                f"received {requested_count}"
            )
        initial = _moments_array(moments)
        configurations = {}
        for group, variants in operations:
            for local_name, transform in variants:
                name = f"{group}/{local_name}"
                try:
                    configurations[name] = _moments_array(transform(initial.copy()))
                except ValueError as error:
                    raise ValueError(f"{name}: {error}") from error
        return configurations

    return count, provider
