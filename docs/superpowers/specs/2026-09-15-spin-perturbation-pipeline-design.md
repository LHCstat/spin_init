# Ordered Spin Perturbation Pipeline Design

## Context

`spin_init` Stage 4 currently creates one unmodified magnetic baseline and
Canting variants for every POSCAR snapshot in `02.disp`. This design extends
that stage with Rotation, Rota_Cant, Random, and Scale operations and allows
users to compose operations in a deterministic order.

The extension remains confined to magnetic-moment generation. Existing
snapshot discovery, `03.spin` directory creation, relative POSCAR/POTCAR
symbolic links, INCAR writing, `MAGMOM`/`M_CONSTR` synchronization,
dpdispatcher submission, and VASP result handling remain unchanged.

## Goals

- Treat `pert_spin` as an ordered sequence of magnetic operations.
- Apply operations in input order and expand every operation over all current
  branches.
- Support Rotation, Canting, Rota_Cant, Random, and Scale.
- Preserve the unmodified `000000` baseline alongside the final leaves.
- Provide stable, descriptive task names that encode the full operation chain.
- Make stochastic operations reproducible through optional non-negative seeds.
- Remain compatible with Python 3.9 and the existing Stage 4 task interface.

## Non-goals

- Changing structural perturbation, AIMD, XDATCAR collection, or VASP
  submission behavior.
- Retaining intermediate operation results as Stage 4 tasks.
- Implementing spin RMSE filtering or DeepMD spin dataset conversion.
- Changing `dpgen init_bulk`, the machine schema, or `pyproject.toml`.

## Input Model

`pert_spin` is an ordered list. Each list item must be a nonempty object with
exactly one operation key. Repeated operation types are allowed because order
is represented by the list rather than JSON object-key order.

```json
"pert_spin": [
  {
    "Rotation": {
      "angle": [45, 90],
      "axis": [[0, 0, 1], [0, 1, 0]]
    }
  },
  {
    "Canting": {
      "angle": [30, 60],
      "seed": 12345
    }
  },
  {
    "Scale": {
      "pert": 0.25,
      "pert_step": 0.05
    }
  }
]
```

An empty or omitted `pert_spin` produces only `000000`.

## Branch Expansion

The pipeline begins with one in-memory branch containing the input magnetic
moments and an empty name. Each operation expands every current branch over
that operation's local variants. Only the final leaves are returned to Stage 4.
The unchanged input is independently retained as `000000` and does not enter
the operation pipeline as an output branch.

For the example above:

- Rotation expands to `2 angles * 2 axes = 4` branches.
- Canting expands every branch over 2 angles, producing 8 branches.
- Scale has 10 nonzero relative deltas, producing 80 final branches.

Expansion order follows input list order, then parameter list order. Cartesian
products use the field order documented for each operation below. This stable
ordering is part of task naming and seeded reproducibility.

## Operation Semantics

All operations accept and return an `N x 3` finite Cartesian moment array in
POSCAR atom order. Unless Scale changes it intentionally, moment magnitude is
preserved. Exact zero moments remain exact zero moments.

### Rotation

```json
{
  "Rotation": {
    "angle": [45, 90],
    "axis": [[0, 0, 1], [0, 1, 0]]
  }
}
```

- `angle` is one finite number or a nonempty list in `[0, 360]` degrees.
- `axis` is one finite, nonzero 3-vector or a nonempty list of such vectors.
- Variants are the Cartesian product `angle * axis`, with angle as the outer
  loop and axis as the inner loop.
- Every nonzero moment is rotated around the normalized global axis with
  Rodrigues' rotation formula and the right-hand rule. Looking along the
  positive axis toward the origin, positive angles appear counterclockwise.
- Zero moments remain zero. Zero axes are rejected.

### Canting

```json
{
  "Canting": {
    "angle": [30, 60],
    "seed": 12345
  }
}
```

- `angle` is one finite number or a nonempty list in `[0, 180]` degrees.
- Each angle creates one local variant.
- For every nonzero moment `m`, construct a stable orthonormal basis `e1`,
  `e2` in the plane perpendicular to `m`, independently sample
  `phi ~ Uniform[0, 2*pi)`, and calculate:

  ```text
  m' = |m| [cos(angle) m_hat + sin(angle) (cos(phi) e1 + sin(phi) e2)]
  ```

- `angle=0` returns the input exactly and `angle=180` reverses every nonzero
  moment exactly. These endpoints do not consume random numbers.
- `seed` is optional. When present, it must be a non-negative integer; explicit
  JSON `null` is rejected.

### Rota_Cant

```json
{
  "Rota_Cant": {
    "R_angle": [30, 60],
    "axis": [[0, 1, 0], [0, 0, 1]],
    "C_angle": [45, 90],
    "seed": 12345
  }
}
```

- `R_angle` and `axis` use Rotation validation and semantics.
- `C_angle` and `seed` use Canting validation and semantics.
- Variants are the Cartesian product `R_angle * axis * C_angle` in that order.
- Each variant first applies Rotation, then applies Canting relative to the
  rotated moments.

### Random

```json
{
  "Random": {
    "num": 5,
    "seed": 12345
  }
}
```

- `num` is a positive integer and creates that many local variants.
- For each variant and each nonzero atom, independently sample a direction
  uniformly on the unit sphere and multiply it by that atom's input magnitude.
- Uniform sphere sampling uses `z ~ Uniform[-1, 1]` and
  `phi ~ Uniform[0, 2*pi)`.
- Zero moments remain zero.
- `seed` follows Canting's optional non-negative integer rule.

### Scale

```json
{
  "Scale": {
    "pert": 0.25,
    "pert_step": 0.05
  }
}
```

- `pert` is finite and satisfies `0 < pert < 1`.
- `pert_step` is finite and positive.
- `pert / pert_step` must be an integer within a floating-point tolerance.
- Relative deltas are generated in stable order:

  ```text
  -pert, -pert + pert_step, ..., -pert_step,
  +pert_step, ..., +pert - pert_step, +pert
  ```

- The zero delta is excluded.
- Each variant applies the same delta to all moments:

  ```text
  m' = (1 + delta) m
  ```

- With `pert < 1`, every scale factor remains positive, so nonzero moment
  directions are preserved. Zero moments remain zero.

## Random Number Semantics

Each Canting, Rota_Cant, or Random operation instance owns one random generator.
The generator advances through current parent branches in stable branch order,
then local variants in parameter order, then atoms in POSCAR order. It is not
reseeded for each atom, variant, or snapshot.

Given the same seed, parameters, snapshots, and task order, the complete output
sequence is reproducible. Omitting a seed initializes a new nondeterministic
sequence. Different seeded operation items own independent generators, even if
they use the same numerical seed.

## Task Naming

The unchanged baseline remains `000000`. Local variant prefixes are:

- Rotation: `R1`, `R2`, ...
- Canting: `C1`, `C2`, ...
- Rota_Cant: `RC1`, `RC2`, ...
- Random: `Rand1`, `Rand2`, ...
- Scale: `S1`, `S2`, ...

For composed operations, local names are joined in execution order with `-`,
for example `R1-C2-S4`. Repeated modes are legal, such as `R1-C1-R2`.
Names remain compatible with the existing Stage 4 safe-name validation.

For every `02.disp/.../POSCAR`, Stage 4 therefore writes:

```text
03.spin/.../000000/INCAR
03.spin/.../R1-C1-S1/INCAR
03.spin/.../R1-C1-S2/INCAR
...
```

Every generated INCAR writes identical final vectors to `MAGMOM` and
`M_CONSTR` through the existing INCAR writer.

## Architecture

`dpgen/data/spin_perturb.py` will contain:

- shared finite-array, scalar-list, axis-list, seed, and scale-grid validators;
- pure vector transforms for Rotation, Canting, Random, and Scale;
- a Rota_Cant transform composed from Rotation followed by Canting;
- operation parsers that produce ordered local variants;
- an ordered pipeline expander that returns the existing provider contract:
  a mapping of safe task names to `N x 3` arrays.

`dpgen/data/spin_tasks.py` remains responsible only for snapshot/INCAR input,
baseline insertion, task validation, and filesystem generation. It will not
contain perturbation mathematics.

The internal `spin_pert_numb` compatibility field remains hidden from users and
is computed from the final branch count before Stage 4 task planning.

## Validation and Errors

Validation occurs before Stage 1 or task-directory creation when Stage 4 make
is requested. Errors identify the operation index and mode, plus the invalid
field or atom where applicable. The implementation rejects:

- a `pert_spin` entry with zero or multiple operation keys;
- unknown operation names or parameters;
- empty/non-finite angle lists and out-of-range angles;
- malformed, non-finite, or zero axes;
- boolean, negative, fractional, or null seeds;
- non-positive/non-integer Random `num`;
- invalid Scale ranges, steps, or non-integral `pert / pert_step` ratios;
- non-finite or malformed magnetic-moment arrays;
- branch names or provider counts inconsistent with Stage 4 expectations.

`spin_action=run` continues to validate the parameter structure but uses the
saved `tasks.json`; it does not regenerate perturbations.

## Testing Strategy

Python 3.9 unit tests will cover:

- Rodrigues rotations around normalized axes, including 0/180/360 degrees;
- Rotation angle-axis Cartesian expansion and zero-axis rejection;
- existing Canting angle, magnitude, endpoint, atomwise sampling, and seed
  contracts;
- Rota_Cant order using a noncommutative example and its Cartesian product;
- Random sphere direction, magnitude preservation, zero moments, count, and
  seeded sequence reproduction;
- Scale delta ordering, zero exclusion, relative factors, divisibility, and
  validation;
- sequential pipeline order, repeated operations, final-only behavior, branch
  counts, and exact composite names;
- parameter validation before filesystem mutation;
- integration with the existing Stage 4 provider and INCAR task generation;
- the existing `spin_init`, CLI, XDATCAR, and `init_bulk` regression tests.

Windows will run all pure-Python tests. Real symbolic-link regression tests and
any VASP integration remain Linux-only according to the existing environment
policy.

## Compatibility

The current one-item Canting input remains valid:

```json
"pert_spin": [{"Canting": {"angle": [30, 60], "seed": 12345}}]
```

Multiple Canting blocks previously produced sibling configurations. They will
now intentionally compose in list order because `pert_spin` is formally an
operation pipeline. The already-removed `Rcut` and `direction` parameters remain
invalid. No machine-file, scheduler, CLI, directory-stage, or packaging schema
changes are required.
