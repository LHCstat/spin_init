# Ordered Spin Perturbation Pipeline Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Extend `spin_init` Stage 4 with ordered Rotation, Canting, Rota_Cant, Random, and Scale operations that Cartesian-expand into reproducible, descriptively named final magnetic configurations.

**Architecture:** Keep magnetic mathematics, parameter validation, operation parsing, and branch expansion in `dpgen/data/spin_perturb.py`, behind the existing `(count, provider)` interface. Keep `dpgen/data/spin_tasks.py` responsible for snapshot validation, the `000000` baseline, INCAR generation, symbolic links, and task directories; `dpgen/data/spin_init.py` continues to validate Stage 4 before any earlier stage mutates the filesystem.

**Tech Stack:** Python 3.9, NumPy, `unittest`, pymatgen-backed existing Stage 4 task handling, dpdispatcher-backed existing submission path.

**Spec:** `docs/superpowers/specs/2026-09-15-spin-perturbation-pipeline-design.md`

## Global Constraints

- Python 3.9 compatibility is mandatory; do not use syntax or standard-library APIs introduced after Python 3.9.
- Do not change `pyproject.toml`, the machine-file schema, scheduler behavior, CLI name, or the Stage 1–3 workflow.
- Keep `000000` as the unmodified baseline and emit only final pipeline leaves as perturbed tasks.
- Keep POSCAR and POTCAR as real relative symbolic links; Windows permission failures must not weaken this requirement.
- Preserve the existing provider signature: `provider(moments, count) -> dict[str, numpy.ndarray]`.
- Every `pert_spin` item contains exactly one operation; execute items in list order and allow repeated operation types.
- Use the exact public mode name `Rota_Cant`.
- Do not retain intermediate operation branches.
- Do not add spin RMSE filtering, DeepMD conversion, or any `04.*` stage.
- Windows runs pure-Python tests; Linux runs real symbolic-link regressions and later VASP integration.

## File Map

- Modify `dpgen/data/spin_perturb.py`: validators, pure transforms, operation compilation, ordered branch expansion, provider construction.
- Modify `dpgen/data/arginfo.py`: update only the `pert_spin` documentation string so normalized parameters describe all supported modes.
- Modify `tests/data/test_spin_perturb.py`: mathematical, validation, naming, ordering, expansion, and seeded-reproducibility tests.
- Modify `tests/data/test_spin_tasks.py`: Stage 4 wiring, pre-filesystem validation, final task-name integration, and run-only validation tests.
- Modify `tests/data/test_spin_init.py`: normalized parameter coverage for the ordered list schema without adding dargs sub-schema restrictions.
- Modify `examples/init/spin-init.json`: one valid ordered pipeline example.
- Modify `SPIN_INIT_GUIDE.md`, `doc/init/spin-init.md`, and `doc/init/spin-init-usage.md`: input schema, formulas, naming, expansion, seeds, and output examples.
- Do not modify `dpgen/data/spin_tasks.py` or `dpgen/data/spin_init.py` unless a failing integration test exposes an incompatibility with the already-defined provider contract.

---

### Task 1: Rotation Transform and Rotation Variants

**Files:**
- Modify: `dpgen/data/spin_perturb.py`
- Test: `tests/data/test_spin_perturb.py`

**Interfaces:**
- Consumes: `_moments_array(moments) -> numpy.ndarray`.
- Produces: `rotate_moments(moments, angle, axis) -> numpy.ndarray`; internal `_rotation_variants(parameters, location) -> list[tuple[str, callable]]` with names `R1`, `R2`, ... in angle-outer/axis-inner order.

- [ ] **Step 1: Add failing Rotation mathematics tests**

```python
class TestRotation(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb
        return spin_perturb

    def test_right_hand_rotation_uses_normalized_global_axis(self):
        moments = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        actual = self.module().rotate_moments(moments, 90, [0, 0, 2])
        np.testing.assert_allclose(actual, [[0, 1, 0], [0, 0, 0]], atol=1e-14)

    def test_rotation_endpoints_preserve_expected_vectors(self):
        moments = np.array([[1.0, 2.0, 3.0]])
        module = self.module()
        np.testing.assert_array_equal(module.rotate_moments(moments, 0, [0, 0, 1]), moments)
        np.testing.assert_allclose(module.rotate_moments(moments, 360, [0, 0, 1]), moments, atol=1e-14)
        np.testing.assert_allclose(module.rotate_moments([[1, 0, 0]], 180, [0, 0, 1]), [[-1, 0, 0]], atol=1e-14)
```

- [ ] **Step 2: Run the new class and confirm the intended failure**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestRotation -v`

Expected: FAIL because `rotate_moments` is not defined.

- [ ] **Step 3: Implement Rodrigues rotation with shared validators**

```python
def rotate_moments(moments, angle, axis):
    result = _moments_array(moments)
    degrees = _numeric_values(angle, "Rotation angle", 0.0, 360.0)[0]
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
```

Implement `_numeric_values` so booleans, empty lists, malformed dimensions, non-finite values, and values outside the requested inclusive range fail with field-specific messages. Implement `_axis_values` so one 3-vector and a nonempty list of 3-vectors are accepted, normalized, and finite nonzero axes are required.

- [ ] **Step 4: Add failing Cartesian-expansion and validation tests**

```python
def test_rotation_cartesian_product_is_angle_outer_axis_inner(self):
    count, provider = self.module().build_spin_perturbation([
        {"Rotation": {"angle": [90, 180], "axis": [[0, 0, 1], [0, 1, 0]]}}
    ])
    configs = provider(np.array([[1.0, 0.0, 0.0]]), count)
    self.assertEqual(count, 4)
    self.assertEqual(list(configs), ["R1", "R2", "R3", "R4"])
    np.testing.assert_allclose(configs["R1"], [[0, 1, 0]], atol=1e-14)
    np.testing.assert_allclose(configs["R2"], [[0, 0, -1]], atol=1e-14)

def test_rotation_rejects_zero_axis_and_unknown_parameters(self):
    for block, message in [
        ({"Rotation": {"angle": 30, "axis": [0, 0, 0]}}, "axis"),
        ({"Rotation": {"angle": 30, "axis": [0, 0, 1], "seed": 1}}, "seed"),
    ]:
        with self.subTest(block=block), self.assertRaisesRegex(ValueError, message):
            self.module().build_spin_perturbation([block])
```

- [ ] **Step 5: Compile Rotation variants and enable the mode in the builder**

Use `itertools.product(angles, axes)` so variant order is exact. Capture immutable angle/axis values in each transform closure and return one local transform per pair. Require each `pert_spin` block to have exactly one key before dispatching to `_rotation_variants`.

- [ ] **Step 6: Run focused and existing Canting tests**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestRotation tests.data.test_spin_perturb.TestCanting -v`

Expected: PASS, except the old test that lists Rotation as unsupported must first remove that now-valid case.

- [ ] **Step 7: Commit Rotation**

```bash
git add dpgen/data/spin_perturb.py tests/data/test_spin_perturb.py
git commit -m "feat: add Rotation spin perturbations"
```

---

### Task 2: Random Direction Transform and Seeded Variants

**Files:**
- Modify: `dpgen/data/spin_perturb.py`
- Test: `tests/data/test_spin_perturb.py`

**Interfaces:**
- Consumes: `_moments_array(moments) -> numpy.ndarray`, `_seed_value(value, location) -> int`.
- Produces: `randomize_moments(moments, rng=None) -> numpy.ndarray`; internal `_random_variants(parameters, location) -> list[tuple[str, callable]]` named `Rand1`, ... with one RNG owned by the operation instance.

- [ ] **Step 1: Add deterministic low-level Random tests**

Extend `SequenceRng.uniform` to support the scalar calls made by both Canting and Random, then add:

```python
class TestRandom(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb
        return spin_perturb

    def test_uniform_sphere_formula_preserves_each_magnitude_and_zero(self):
        rng = SequenceRng([0.0, 0.0, 0.0, np.pi / 2])
        actual = self.module().randomize_moments(
            [[2, 0, 0], [0, 0, 0], [0, 0, 3]], rng=rng
        )
        np.testing.assert_allclose(actual, [[2, 0, 0], [0, 0, 0], [0, 3, 0]], atol=1e-14)
```

The sequence represents `(z=0, phi=0)` for the first nonzero atom and `(z=0, phi=pi/2)` for the second; the zero atom consumes no samples.

- [ ] **Step 2: Run and confirm failure**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestRandom.test_uniform_sphere_formula_preserves_each_magnitude_and_zero -v`

Expected: FAIL because `randomize_moments` is not defined.

- [ ] **Step 3: Implement uniform sphere sampling**

```python
def randomize_moments(moments, rng=None):
    result = _moments_array(moments)
    if rng is None:
        rng = np.random.default_rng()
    _require_uniform_rng(rng, "Random")
    for atom_index, moment in enumerate(result):
        _, magnitude = _unit_and_magnitude(moment)
        if magnitude == 0.0:
            continue
        if not np.isfinite(magnitude):
            raise ValueError(f"atom {atom_index}: magnetic-moment magnitude is too large")
        z = rng.uniform(-1.0, 1.0)
        phi = rng.uniform(0.0, 2.0 * np.pi)
        radius = np.sqrt(max(0.0, 1.0 - z * z))
        result[atom_index] = magnitude * np.array([radius * np.cos(phi), radius * np.sin(phi), z])
    return result
```

- [ ] **Step 4: Add provider count, sequence, and validation tests**

```python
def test_random_num_creates_seeded_sequence_without_reseeding(self):
    params = [{"Random": {"num": 3, "seed": 789}}]
    count_a, provider_a = self.module().build_spin_perturbation(params)
    count_b, provider_b = self.module().build_spin_perturbation(params)
    moments = np.array([[0.0, 0.0, 2.0], [0.0, 0.0, 0.0]])
    first_a = provider_a(moments, count_a)
    next_a = provider_a(moments, count_a)
    first_b = provider_b(moments, count_b)
    self.assertEqual(list(first_a), ["Rand1", "Rand2", "Rand3"])
    for name in first_a:
        np.testing.assert_array_equal(first_a[name], first_b[name])
        np.testing.assert_allclose(np.linalg.norm(first_a[name][0]), 2.0)
        self.assertFalse(np.array_equal(first_a[name], next_a[name]))

def test_random_requires_positive_integer_num_and_valid_seed(self):
    for value in [0, -1, 1.5, True, None]:
        with self.subTest(value=value), self.assertRaises(ValueError):
            self.module().build_spin_perturbation([{"Random": {"num": value}}])
```

- [ ] **Step 5: Compile Random variants using one operation-owned RNG**

Create `rng = np.random.default_rng(seed)` once inside `_random_variants`; every local closure refers to that same RNG. Validate exactly `num` and optional `seed`; reject unknown fields. The provider iterates variants in `Rand1` to `RandN` order and never resets the RNG between provider calls.

- [ ] **Step 6: Run Random and Canting seed regressions**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestRandom tests.data.test_spin_perturb.TestCanting.test_seed_reproduces_the_complete_configuration_sequence -v`

Expected: PASS.

- [ ] **Step 7: Commit Random**

```bash
git add dpgen/data/spin_perturb.py tests/data/test_spin_perturb.py
git commit -m "feat: add Random spin perturbations"
```

---

### Task 3: Relative Scale Grid and Scale Variants

**Files:**
- Modify: `dpgen/data/spin_perturb.py`
- Test: `tests/data/test_spin_perturb.py`

**Interfaces:**
- Consumes: `_moments_array(moments) -> numpy.ndarray`.
- Produces: `scale_moments(moments, delta) -> numpy.ndarray`; internal `_scale_deltas(pert, pert_step) -> list[float]`; internal `_scale_variants(parameters, location) -> list[tuple[str, callable]]` named `S1`, ... in negative-then-positive delta order.

- [ ] **Step 1: Add failing Scale grid and transform tests**

```python
class TestScale(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb
        return spin_perturb

    def test_grid_excludes_zero_and_scale_is_relative(self):
        module = self.module()
        np.testing.assert_allclose(module._scale_deltas(0.25, 0.05), [
            -0.25, -0.20, -0.15, -0.10, -0.05,
             0.05,  0.10,  0.15,  0.20,  0.25,
        ], rtol=0, atol=1e-15)
        np.testing.assert_allclose(module.scale_moments([[2, -4, 0]], -0.25), [[1.5, -3, 0]])

    def test_scale_provider_names_follow_delta_order(self):
        count, provider = self.module().build_spin_perturbation([
            {"Scale": {"pert": 0.10, "pert_step": 0.05}}
        ])
        configs = provider(np.array([[0.0, 0.0, 2.0]]), count)
        self.assertEqual(count, 4)
        self.assertEqual(list(configs), ["S1", "S2", "S3", "S4"])
        np.testing.assert_allclose([configs[name][0, 2] for name in configs], [1.8, 1.9, 2.1, 2.2])
```

- [ ] **Step 2: Run and confirm failure**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestScale -v`

Expected: FAIL because the Scale functions are not defined or the mode is unsupported.

- [ ] **Step 3: Implement a stable integer-indexed grid**

```python
def _scale_deltas(pert, pert_step):
    pert_value = _single_finite_number(pert, "Scale pert")
    step_value = _single_finite_number(pert_step, "Scale pert_step")
    if not 0.0 < pert_value < 1.0:
        raise ValueError("Scale pert must satisfy 0 < pert < 1")
    if step_value <= 0.0:
        raise ValueError("Scale pert_step must be positive")
    steps = int(round(pert_value / step_value))
    if steps < 1 or not np.isclose(steps * step_value, pert_value, rtol=1e-12, atol=1e-15):
        raise ValueError("Scale pert / pert_step must be an integer")
    negative = [-step_value * index for index in range(steps, 0, -1)]
    positive = [step_value * index for index in range(1, steps + 1)]
    return [float(value) for value in negative + positive]

def scale_moments(moments, delta):
    return _moments_array(moments) * (1.0 + float(delta))
```

Generate values by integer multiplication, not `numpy.arange`, to prevent an accidental zero or missing endpoint.

- [ ] **Step 4: Add explicit validation tests**

```python
def test_scale_rejects_invalid_range_step_ratio_and_fields(self):
    invalid = [
        {"pert": 0}, {"pert": 1}, {"pert": -0.1},
        {"pert": 0.25, "pert_step": 0},
        {"pert": 0.25, "pert_step": 0.06},
        {"pert": 0.25, "pert_step": 0.05, "seed": 1},
    ]
    for parameters in invalid:
        with self.subTest(parameters=parameters), self.assertRaises(ValueError):
            self.module().build_spin_perturbation([{"Scale": parameters}])
```

- [ ] **Step 5: Compile Scale variants and preserve zero vectors**

Validate exactly `pert` and `pert_step`, call `_scale_deltas`, and capture each delta in its `S#` transform. `_moments_array` plus multiplication keeps exact zero moments unchanged and rejects malformed/non-finite input arrays.

- [ ] **Step 6: Run Scale tests**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestScale -v`

Expected: PASS.

- [ ] **Step 7: Commit Scale**

```bash
git add dpgen/data/spin_perturb.py tests/data/test_spin_perturb.py
git commit -m "feat: add relative Scale spin perturbations"
```

---

### Task 4: Rota_Cant Composite Operation

**Files:**
- Modify: `dpgen/data/spin_perturb.py`
- Test: `tests/data/test_spin_perturb.py`

**Interfaces:**
- Consumes: `rotate_moments(moments, angle, axis)`, `cant_moments(moments, angle, rng)`, Rotation/Canting validators.
- Produces: `rota_cant_moments(moments, rotation_angle, axis, canting_angle, rng=None) -> numpy.ndarray`; internal `_rota_cant_variants(parameters, location) -> list[tuple[str, callable]]` named `RC1`, ... in `R_angle`-outer/`axis`-middle/`C_angle`-inner order with Rotation always applied before Canting.

- [ ] **Step 1: Add a failing noncommutative order test**

```python
class TestRotaCant(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb
        return spin_perturb

    def test_rotation_is_applied_before_canting(self):
        rotated = self.module().rotate_moments([[1.0, 0.0, 0.0]], 90, [0, 0, 1])
        expected = self.module().cant_moments(rotated, 90, rng=SequenceRng([0.0]))
        actual = self.module().rota_cant_moments(
            [[1.0, 0.0, 0.0]],
            rotation_angle=90,
            axis=[0, 0, 1],
            canting_angle=90,
            rng=SequenceRng([0.0]),
        )
        np.testing.assert_allclose(actual, expected, atol=1e-14)
```

- [ ] **Step 2: Run and confirm failure**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestRotaCant.test_rotation_is_applied_before_canting -v`

Expected: FAIL because `Rota_Cant` is unsupported.

- [ ] **Step 3: Implement composition by reusing the two transforms**

```python
def rota_cant_moments(moments, rotation_angle, axis, canting_angle, rng=None):
    rotated = rotate_moments(moments, rotation_angle, axis)
    return cant_moments(rotated, canting_angle, rng=rng)
```

Validate exactly `R_angle`, `axis`, `C_angle`, and optional `seed`. Create one RNG for the entire Rota_Cant operation instance, then build variants with `itertools.product(rotation_angles, axes, canting_angles)`.

- [ ] **Step 4: Add Cartesian count, naming, and seed tests**

```python
def test_rota_cant_cartesian_product_and_seed_are_stable(self):
    parameters = [{"Rota_Cant": {
        "R_angle": [30, 60],
        "axis": [[0, 1, 0], [0, 0, 1]],
        "C_angle": [45, 90],
        "seed": 12345,
    }}]
    count_a, provider_a = self.module().build_spin_perturbation(parameters)
    count_b, provider_b = self.module().build_spin_perturbation(parameters)
    moments = np.array([[0.0, 0.0, 2.0]])
    configs_a = provider_a(moments, count_a)
    configs_b = provider_b(moments, count_b)
    self.assertEqual(count_a, 8)
    self.assertEqual(list(configs_a), [f"RC{i}" for i in range(1, 9)])
    for name in configs_a:
        np.testing.assert_array_equal(configs_a[name], configs_b[name])
        np.testing.assert_allclose(np.linalg.norm(configs_a[name][0]), 2.0)
```

- [ ] **Step 5: Add field-specific invalid-input cases**

Test missing `R_angle`, missing `axis`, missing `C_angle`, zero axis, out-of-range Rotation/Canting angles, invalid/null seed, and unknown fields. Assert the error includes `pert_spin[0].Rota_Cant` plus the invalid field.

- [ ] **Step 6: Run Rota_Cant and component regressions**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestRotaCant tests.data.test_spin_perturb.TestRotation tests.data.test_spin_perturb.TestCanting -v`

Expected: PASS.

- [ ] **Step 7: Commit Rota_Cant**

```bash
git add dpgen/data/spin_perturb.py tests/data/test_spin_perturb.py
git commit -m "feat: add Rota_Cant spin perturbations"
```

---

### Task 5: Ordered Pipeline Expansion, Repeated Modes, and Composite Names

**Files:**
- Modify: `dpgen/data/spin_perturb.py`
- Test: `tests/data/test_spin_perturb.py`

**Interfaces:**
- Consumes: each mode compiler returning local `(name, transform)` variants.
- Produces: `build_spin_perturbation(pert_spin) -> tuple[int, Optional[callable]]`, where the provider returns only final leaves in stable order and keeps operation-owned RNG state across snapshots.

- [ ] **Step 1: Replace the old sibling-Canting expectation with ordered composition tests**

```python
class TestOrderedPipeline(unittest.TestCase):
    def module(self):
        from dpgen.data import spin_perturb
        return spin_perturb

    def test_operations_expand_current_branches_and_emit_only_final_leaves(self):
        count, provider = self.module().build_spin_perturbation([
            {"Rotation": {"angle": [90, 180], "axis": [0, 0, 1]}},
            {"Scale": {"pert": 0.1, "pert_step": 0.1}},
        ])
        configs = provider(np.array([[1.0, 0.0, 0.0]]), count)
        self.assertEqual(count, 4)
        self.assertEqual(list(configs), ["R1-S1", "R1-S2", "R2-S1", "R2-S2"])
        self.assertNotIn("R1", configs)
        np.testing.assert_allclose(configs["R1-S1"], [[0, 0.9, 0]], atol=1e-14)
        np.testing.assert_allclose(configs["R1-S2"], [[0, 1.1, 0]], atol=1e-14)

    def test_repeated_modes_compose_in_list_order(self):
        count, provider = self.module().build_spin_perturbation([
            {"Rotation": {"angle": 90, "axis": [0, 0, 1]}},
            {"Canting": {"angle": 0, "seed": 1}},
            {"Rotation": {"angle": 90, "axis": [0, 0, 1]}},
        ])
        configs = provider(np.array([[1.0, 0.0, 0.0]]), count)
        self.assertEqual(list(configs), ["R1-C1-R1"])
        np.testing.assert_allclose(configs["R1-C1-R1"], [[-1, 0, 0]], atol=1e-14)
```

Local numbering restarts for each operation instance, so the second Rotation's only local name is `R1`; the complete chain remains unambiguous by position.

- [ ] **Step 2: Run and confirm the pipeline tests fail**

Run: `python3.9 -m unittest tests.data.test_spin_perturb.TestOrderedPipeline -v`

Expected: FAIL because the existing builder flattens blocks as sibling configurations.

- [ ] **Step 3: Implement compile-time operation metadata and provider-time expansion**

```python
def build_spin_perturbation(pert_spin):
    operations = _compile_operations(pert_spin)
    if not operations:
        return 0, None
    count = 1
    for variants in operations:
        count *= len(variants)

    def provider(moments, requested_count):
        if requested_count != count:
            raise ValueError(f"spin provider expected {count} configurations, received {requested_count}")
        branches = [("", _moments_array(moments))]
        for variants in operations:
            expanded = []
            for parent_name, parent_moments in branches:
                for local_name, transform in variants:
                    name = local_name if not parent_name else f"{parent_name}-{local_name}"
                    expanded.append((name, _moments_array(transform(parent_moments))))
            branches = expanded
        return {name: values for name, values in branches}

    return count, provider
```

Keep RNGs inside compiled operation transforms, so each stochastic operation instance owns one generator. Stable iteration order is parent branch, then local variant, then atom.

- [ ] **Step 4: Enforce the exact one-key block and mode dispatch contract**

```python
def test_each_pipeline_item_has_exactly_one_known_operation(self):
    invalid = [
        {},
        {"Rotation": {"angle": 30, "axis": [0, 0, 1]}, "Scale": {"pert": 0.1, "pert_step": 0.1}},
        {"rotation": {"angle": 30, "axis": [0, 0, 1]}},
        "Canting",
    ]
    for index, block in enumerate(invalid):
        with self.subTest(block=block), self.assertRaises((ValueError, NotImplementedError)):
            self.module().build_spin_perturbation([block])
```

Use a dispatch table with the exact keys `Rotation`, `Canting`, `Rota_Cant`, `Random`, and `Scale`. Error text includes `pert_spin[index]` and the mode when available.

- [ ] **Step 5: Test operation-owned RNG order across parent branches and snapshots**

Build identical seeded providers twice for Rotation → Canting and call each provider on two snapshots. Assert corresponding arrays are identical between providers, later snapshots differ from earlier snapshots for non-endpoint Canting, and result key order is identical. Add the equivalent Random check. This verifies no reseed per parent, variant, atom, or snapshot.

- [ ] **Step 6: Run all perturbation tests**

Run: `python3.9 -m unittest tests.data.test_spin_perturb -v`

Expected: PASS. Update `test_multiple_canting_blocks_continue_configuration_numbering` to assert composition names such as `C1-C1`, `C1-C2` and count `2`, rather than the old sibling count `3`.

- [ ] **Step 7: Commit the ordered pipeline**

```bash
git add dpgen/data/spin_perturb.py tests/data/test_spin_perturb.py
git commit -m "feat: compose ordered spin perturbation operations"
```

---

### Task 6: Stage 4 Integration and Pre-Filesystem Validation

**Files:**
- Modify: `tests/data/test_spin_tasks.py`
- Modify: `tests/data/test_spin_init.py`
- Modify only if required by a failing test: `dpgen/data/spin_init.py`

**Interfaces:**
- Consumes: `build_spin_perturbation` provider count and final-name mapping; existing `make_spin_tasks(jdata, mdata, perturb=None)`.
- Produces: Stage 4 `spin_pert_numb` equal to the final branch product, with `000000` added only by `plan_spin_tasks`.

- [ ] **Step 1: Replace unsupported-Rotation tests with ordered-pipeline wiring tests**

```python
def test_stage4_wires_final_pipeline_count_and_names_into_task_generation(self):
    module = self.module()
    self.jdata["spin_action"] = "make"
    self.jdata["pert_spin"] = [
        {"Rotation": {"angle": [90, 180], "axis": [0, 0, 1]}},
        {"Scale": {"pert": 0.1, "pert_step": 0.1}},
    ]
    param = self.root / "input.json"
    param.write_text(json.dumps(self.jdata))
    with mock.patch.object(module, "make_spin_tasks") as make:
        spin_init.gen_spin_init(argparse.Namespace(PARAM=str(param), MACHINE=None))
    generated_jdata, _ = make.call_args.args
    provider = make.call_args.kwargs["perturb"]
    self.assertEqual(generated_jdata["spin_pert_numb"], 4)
    self.assertEqual(list(provider(np.array([[1.0, 0.0, 0.0]]), 4)), [
        "R1-S1", "R1-S2", "R2-S1", "R2-S2"
    ])
```

- [ ] **Step 2: Add a full Stage 4 planning test with the baseline**

Use the existing temporary POSCAR/spin-INCAR fixture, configure a one-variant Rotation followed by two Scale variants, call `plan_spin_tasks`, and assert the task suffixes are exactly `000000`, `R1-S1`, and `R1-S2` for every snapshot. Assert generated moment arrays have the expected values and MAGMOM/M_CONSTR writing remains covered by the existing task tests.

- [ ] **Step 3: Update early-validation tests**

Replace the old valid Rotation failure case with an actually invalid item such as:

```python
self.jdata["pert_spin"] = [{"Rotation": {"angle": 30, "axis": [0, 0, 0]}}]
```

For `stages=[1, 4]`, assert `make_spin_init_structures` is not called and the exception contains `pert_spin[0].Rotation.axis`. For `spin_action=run`, assert the same parameter validation occurs before submission even though perturbations are not regenerated.

- [ ] **Step 4: Add normalized schema coverage**

In `tests/data/test_spin_init.py`, normalize a `pert_spin` list containing all five exact mode keys and assert list order and nested values survive unchanged. This intentionally keeps `Argument("pert_spin", list[dict])` permissive; mode-specific validation remains in `build_spin_perturbation`, where errors can include operation indexes.

- [ ] **Step 5: Run integration tests**

Run: `python3.9 -m unittest tests.data.test_spin_tasks tests.data.test_spin_init -v`

Expected: PASS on Windows except existing real-symlink cases that are explicitly skipped by the test environment. If the existing `spin_init.py` provider wiring already passes, do not edit it.

- [ ] **Step 6: Commit Stage 4 integration tests**

```bash
git add tests/data/test_spin_tasks.py tests/data/test_spin_init.py
git commit -m "test: cover ordered spin pipeline integration"
```

---

### Task 7: User-Facing Parameters, Examples, and Documentation

**Files:**
- Modify: `dpgen/data/arginfo.py`
- Modify: `examples/init/spin-init.json`
- Modify: `SPIN_INIT_GUIDE.md`
- Modify: `doc/init/spin-init.md`
- Modify: `doc/init/spin-init-usage.md`
- Test: `tests/data/test_spin_init.py`

**Interfaces:**
- Consumes: exact mode and field names implemented in Tasks 1–5.
- Produces: copyable JSON examples and documentation matching runtime validation and output naming.

- [ ] **Step 1: Update the concise arginfo description**

Use this exact meaning without introducing a nested dargs schema:

```python
doc=(
    "Ordered magnetic perturbation pipeline. Each item contains exactly one "
    "of Rotation, Canting, Rota_Cant, Random, or Scale."
),
```

- [ ] **Step 2: Replace the example with an ordered pipeline**

Use a small, tractable example rather than the full 80-way product:

```json
"pert_spin": [
  {"Rotation": {"angle": 45, "axis": [0, 0, 1]}},
  {"Canting": {"angle": [30, 60], "seed": 12345}},
  {"Scale": {"pert": 0.10, "pert_step": 0.10}}
]
```

This creates `1 * 2 * 2 = 4` final perturbations plus `000000`.

- [ ] **Step 3: Document the shared ordered expansion contract**

In all three guides, state:

- list order is execution order;
- each item has exactly one mode;
- repeated modes are allowed;
- each operation expands every current branch by Cartesian product;
- only final leaves are written;
- `000000` is always the unchanged baseline;
- composite names follow execution order, such as `R1-C2-S1`.

- [ ] **Step 4: Document every operation with formulas and validation**

Include the exact fields and formulas from the approved spec:

- Rotation: `angle`, `axis`, Rodrigues/right-hand rule, angle range `[0,360]`.
- Canting: `angle`, optional `seed`, atomwise random azimuth, angle range `[0,180]`.
- Rota_Cant: `R_angle`, `axis`, `C_angle`, optional `seed`, Rotation then Canting.
- Random: positive integer `num`, optional `seed`, uniform sphere directions, unchanged norms.
- Scale: `pert`, `pert_step`, relative `(1+delta)m`, `0 < pert < 1`, integral ratio, zero delta excluded.

State that zero moments stay zero and stochastic RNG state advances across variants, parent branches, and snapshots.

- [ ] **Step 5: Document output counts and names with one calculated example**

For Rotation with 2 angles × 2 axes, Canting with 2 angles, and Scale with `pert=0.25`, `pert_step=0.05`, show `4 * 2 * 10 = 80` final perturbed INCAR tasks plus `000000` per snapshot. Show representative paths:

```text
03.spin/.../000000/INCAR
03.spin/.../R1-C1-S1/INCAR
03.spin/.../R4-C2-S10/INCAR
```

- [ ] **Step 6: Validate JSON and documentation references**

Run: `python3.9 -m json.tool examples/init/spin-init.json`

Run: `python3.9 -m unittest tests.data.test_spin_init tests.test_cli -v`

Expected: JSON parses and tests PASS.

- [ ] **Step 7: Commit documentation**

```bash
git add dpgen/data/arginfo.py examples/init/spin-init.json SPIN_INIT_GUIDE.md doc/init/spin-init.md doc/init/spin-init-usage.md tests/data/test_spin_init.py
git commit -m "docs: describe ordered spin perturbation modes"
```

---

### Task 8: Python 3.9 Regression Verification and Linux Handoff

**Files:**
- Verify: all files changed in Tasks 1–7
- Preserve untracked: `tests/data/almg.fcc.01x01x01/`

**Interfaces:**
- Consumes: complete ordered perturbation pipeline.
- Produces: recorded Windows pure-Python evidence and exact Linux verification commands; no behavioral workaround for Windows symlink permissions.

- [ ] **Step 1: Run the focused Python 3.9 suite on Windows**

Run from this checkout using the known Python 3.9 interpreter:

```powershell
C:\Users\Murphy\Documents\ChatGPT\dpgen\dpgen_venv\python.exe -m unittest tests.data.test_spin_perturb tests.data.test_spin_tasks tests.data.test_spin_init tests.test_cli tests.data.test_spin_xdatcar -v
```

Expected: all pure-Python tests PASS; only tests explicitly requiring real Windows symbolic-link privilege may SKIP.

- [ ] **Step 2: Run upstream initialization regressions on Windows where supported**

```powershell
C:\Users\Murphy\Documents\ChatGPT\dpgen\dpgen_venv\python.exe -m unittest tests.data.test_gen_bulk tests.data.test_coll_vasp -v
```

Expected: `test_coll_vasp` and non-symlink `test_gen_bulk` cases PASS; `WinError 1314` remains an environment limitation and must not be fixed by copying links or modifying upstream behavior.

- [ ] **Step 3: Check formatting and repository scope**

Run: `python3.9 -m compileall -q dpgen tests`

Run: `git diff --check`

Run: `git status --short`

Expected: compilation and whitespace checks succeed; only intended tracked files plus the pre-existing untracked `tests/data/almg.fcc.01x01x01/` appear.

- [ ] **Step 4: Run Linux Python 3.9 tests in the supercomputer checkout**

From the Linux checkout under `lhc/dpgenmake`, run:

```bash
python3.9 -m unittest tests.data.test_spin_perturb tests.data.test_spin_tasks tests.data.test_spin_init tests.test_cli tests.data.test_spin_xdatcar tests.data.test_gen_bulk tests.data.test_coll_vasp -v
```

Expected: all tests PASS, including real POSCAR/POTCAR symbolic-link assertions.

- [ ] **Step 5: Defer VASP integration until unit and Linux regression gates pass**

Use the already-defined minimal integration input with `scale=[1.0]`, `pert_numb=1`, and `md_nstep=2` or `3`. Verify separately that MD completion is recognized, `XDATCAR` is backwarded, snapshots populate `02.disp`, and Stage 4 creates/optionally submits final magnetic tasks. Merely placing `vasp.slurm` in a task directory is not evidence that it was submitted with `sbatch`; dpdispatcher submission logs are the submission evidence.

- [ ] **Step 6: Request an independent code review and fix only evidenced issues**

Review against the approved spec, emphasizing operation order, Cartesian counts, seed lifetime, exact names, Python 3.9 syntax, and unchanged Stage 1–3 behavior. Re-run the focused suite after every review fix.

- [ ] **Step 7: Commit any verification-only corrections, then push main**

If verification required tracked corrections:

```bash
git add dpgen/data/spin_perturb.py dpgen/data/arginfo.py tests/data/test_spin_perturb.py tests/data/test_spin_tasks.py tests/data/test_spin_init.py examples/init/spin-init.json SPIN_INIT_GUIDE.md doc/init/spin-init.md doc/init/spin-init-usage.md
git commit -m "fix: address ordered spin pipeline verification"
```

Then confirm `git status --short`, preserve the generated untracked test directory, and push the `main` branch only after all applicable gates pass.
