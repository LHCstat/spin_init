# `spin_init` Detailed Usage Reference

English | [Chinese](spin-init-usage.md)

## Purpose and workflow

`spin_init` starts from an existing POSCAR, creates structurally perturbed VASP
AIMD tasks, splits every XDATCAR trajectory into independent structures, and
builds noncollinear magnetic static or relaxation calculations for each
snapshot.

```text
POSCAR
  -> supercell, scaling, and box/atomic perturbation
  -> VASP AIMD
  -> XDATCAR -> independent POSCAR snapshots
  -> noncollinear magnetic static or relaxation VASP tasks
  -> completed magnetic calculations -> DeepMD magnetic data
```

Run it with:

```bash
dpgen spin_init PARAM [MACHINE]
```

The five stages are:

| Stage | Purpose | Output directory |
| --- | --- | --- |
| 1 | Build the supercell, scale it, and generate structural perturbations. | `00.scale_pert` |
| 2 | Create AIMD tasks and submit them when MACHINE is provided. | `01.md` |
| 3 | Check AIMD outputs, parse XDATCAR, and export snapshots. | `02.disp` |
| 4 | Create and/or submit magnetic static or relaxation calculations for every snapshot. | `03.spin` |
| 5 | Validate and collect all magnetic results and export DeepMD data. | `04.data` |

## Input files

The complete workflow uses two INCAR roles. Stage 4 accepts either one magnetic
template or an ordered list of templates.

| File | Purpose |
| --- | --- |
| `POSCAR` | Initial crystal structure. |
| `INCAR.md` | AIMD settings for Stage 2. |
| `INCAR.spin` or several `INCAR.state_*` files | Noncollinear static or relaxation templates for Stage 4. |
| `POTCAR` or several POTCAR fragments | VASP pseudopotentials in POSCAR element order. |
| `spin-init.json` | Workflow parameters. |
| `machine.json` | DPDispatcher configuration for automatic submission. |
| `KPOINTS`, etc. | Additional inputs supplied through `user_forward_files`. |

### Writing `INCAR.md`

This file controls AIMD. A minimal short test may use:

```text
SYSTEM = spin_init_md
ENCUT = 500
EDIFF = 1E-5
IBRION = 0
NSW = 3
POTIM = 1.0
TEBEG = 300
TEEND = 300
SMASS = 0
ISIF = 2
ISMEAR = 1
SIGMA = 0.1
LWAVE = F
LCHARG = F
```

When `NSW` differs from `md_nstep` in PARAM, the workflow uses `NSW` as the
expected AIMD step count.

### Writing `INCAR.spin`

This is an independent Stage 4 template; it neither replaces nor modifies
`INCAR.md`. Every template must satisfy the following rules:

- both `MAGMOM` and `M_CONSTR` are present and initially identical;
- both fields contain exactly three finite Cartesian components per POSCAR
  atom;
- noncollinear magnetism is enabled with `LNONCOLLINEAR = .TRUE.` or `LSORBIT`;
- `NSW >= 0`; user values for `NSW`, `IBRION`, `ISIF`, and other VASP settings
  are retained.

A two-atom static example is:

```text
SYSTEM = spin_init_static
ENCUT = 500
EDIFF = 1E-6
NSW = 0
IBRION = -1
LNONCOLLINEAR = .TRUE.
I_CONSTRAINED_M = 2
LAMBDA = 10
MAGMOM = 0 0 2  0 0 0
M_CONSTR = 0 0 2  0 0 0
LWAVE = .FALSE.
LCHARG = .FALSE.
```

VASP repetition syntax and backslash continuations, such as `6*0.0`, are
supported. The workflow parses and validates the moments, writes a private
INCAR for every task, and keeps that task's `M_CONSTR` synchronized with its
perturbed `MAGMOM`.

### Stage 4 structural and cell relaxation

To relax a structure, replace the static settings with values such as:

```text
NSW = 50
IBRION = 2
ISIF = 3
EDIFFG = -0.02
```

`NSW` is the maximum number of ionic steps; a converged run may stop earlier.
`IBRION=1`, `2`, or `3` with `NSW>0` denotes common relaxation algorithms.
`ISIF=2` relaxes positions at a fixed cell, while `ISIF=3` permits changes to
positions, cell shape, and volume. These values are examples and must be
validated for the physical system. If `IBRION` is omitted while `NSW>0`, VASP
uses its own default rather than having the workflow inject `IBRION=2`.

See the VASP documentation for
[IBRION](https://vasp.at/wiki/index.php/IBRION),
[ISIF](https://vasp.at/wiki/index.php/ISIF), and
[NSW](https://vasp.at/wiki/index.php/NSW).

Relaxation submissions automatically request `CONTCAR`. Static and relaxation
templates may be mixed in one submission. Because DPDispatcher uses one shared
backward-file list, the presence of any relaxation task requests `CONTCAR` for
the entire submission, while final-structure validation is applied only to
relaxation tasks. User backward files are appended without duplicates.

Result validation separates two questions:

- **Normal termination:** OUTCAR contains a timing footer and force output, and
  OSZICAR is nonempty. A static `NSW=0` task retains the single-force-block
  check; an ionic calculation may contain several blocks and may stop before
  its maximum `NSW`.
- **Structural convergence:** a relaxation task also needs a nonempty,
  parseable CONTCAR with the expected atom count and species order. Lattice and
  coordinates may change. VASP's `reached required accuracy - stopping
  structural energy minimisation` marker confirms structural convergence. A
  normally terminated relaxation without that marker produces an explicit
  warning rather than a convergence claim.

`check_spin_results` returns the number of normally terminated tasks, including
relaxations for which structural convergence was not confirmed. This count is
not an electronic or magnetic convergence statement. CONTCAR represents the
last ionic structure and can therefore be unconverged. The initial POSCAR
remains a relative symlink to `02.disp`; neither that link nor its snapshot
source is replaced. See the VASP
[CONTCAR documentation](https://vasp.at/wiki/index.php/CONTCAR).

## Writing `spin-init.json`

This example runs Stages 1 through 4. Append `5` to `stages` when data export is
required.

```json
{
  "stages": [1, 2, 3, 4],
  "from_poscar_path": "./POSCAR",
  "out_dir": "./run_spin",
  "super_cell": [1, 1, 1],
  "scale": [1.0],
  "pert_numb": 1,
  "pert_box": 0.03,
  "pert_atom": 0.01,
  "md_incar": "./INCAR.md",
  "md_nstep": 3,
  "spin_incar": "./INCAR.spin",
  "pert_spin": [
    {"Rotation": {"angle": 45, "axis": [0, 0, 1]}},
    {"Canting": {"angle": [30, 60], "seed": 12345}},
    {"Scale": {"pert": 0.10, "pert_step": 0.10}}
  ],
  "spin_action": "make_run",
  "potcars": ["./POTCAR"]
}
```

Stage 4 fields are:

| Field | Meaning |
| --- | --- |
| `spin_incar` | One magnetic static/relaxation INCAR path, or a nonempty ordered list of paths. |
| `pert_spin` | Independent operation groups in input order. Supported modes are `Rotation`, `Canting`, `Rota_Cant`, `Random`, and `Scale`. |
| `spin_pert_numb` | Internal compatibility field; omit it or leave it at `0`. |
| `spin_action` | `make`, `run`, or `make_run`. |

`spin_action="make"` creates `03.spin` without submission. `run` submits an
existing `03.spin` and requires MACHINE. `make_run` creates and submits tasks
when MACHINE is present, or only creates them when MACHINE is omitted.

### One or multiple initial magnetic states

The string form uses one template and does not add an INCAR index layer:

```json
"spin_incar": "./INCAR.spin"
```

To apply the same operations to multiple initial states, use an ordered list:

```json
"spin_incar": [
  "./INCAR.state_1",
  "./INCAR.state_2"
]
```

The list must be nonempty, every entry must be a nonempty string, and resolved
paths must be unique. List input adds `000-incar`, `001-incar`, ... between the
snapshot and its magnetic groups. A one-item list still adds `000-incar`.

Each template supplies its own initial `MAGMOM/M_CONSTR`, after which the same
`pert_spin` list is applied. Traversal order is snapshot, INCAR input order,
operation input order, and local variant order. All templates share each
stochastic operation's continuous random-number sequence. Their random results
therefore differ, while identical seeds and input order reproduce the complete
sequence.

Every snapshot/template pair has one unchanged baseline at
`origin/000000`. `origin` is a category alongside operation groups such as
`000-rotation`. Each `pert_spin` item must contain exactly one mode; repeated
modes are allowed. Different items independently operate on the original
template moments. Cartesian products occur only among parameter lists inside
the same item, so variant counts from different items are added rather than
multiplied.

Groups are named `<zero-based-input-index>-<lowercase-mode>`, for example
`000-rotation`, `001-canting`, and `002-scale`. Local names are `R#`, `C#`,
`RC#`, `Rand#`, and `S#`. Repeating a mode creates another indexed group rather
than composing with the earlier one. Use `Rota_Cant` for the explicit sequence
of Rotation followed by Canting.

### Canting parameters and definition

`angle` is the angle in degrees between the input and perturbed moments, in the
closed interval `[0, 180]`. It may be one number or a list; `[30, 60]` generates
`C1` and `C2`.

For every nonzero input moment **a**, the workflow constructs an orthonormal
basis **e1**, **e2** in the plane normal to **a** and independently samples an
azimuth `phi` for each atom from `[0, 2*pi)`:

```text
a' = |a| [cos(angle) a_hat + sin(angle) (cos(phi) e1 + sin(phi) e2)]
```

Thus the angle between **a'** and **a** is exactly `angle`, and the magnitude
remains `|a|`. Zero moments remain zero. `angle=0` preserves the vector exactly;
`angle=180` reverses each nonzero vector exactly. Those endpoints consume no
random values.

`seed` is optional and, when present, must be a non-negative integer. The same
seed, inputs, and task order reproduce the complete sequence. An omitted seed
produces a new random sequence. `Rcut` and `direction` are not Canting fields.

### Rotation, Rota_Cant, Random, and Scale

`Rotation.angle` is one value or a nonempty list in `[0, 360]`. `axis` is one
nonzero three-vector or a list of vectors. Local variants follow
`angle x axis` order. Rodrigues' formula applies a right-hand-rule rotation
around the normalized global axis. Zero vectors remain zero and nonzero
magnitudes are preserved.

```json
{"Rotation": {"angle": [45, 90], "axis": [[0, 0, 1], [0, 1, 0]]}}
```

`Rota_Cant` always applies Rotation first and Canting second relative to the
rotated moments. It generates `R_angle x axis x C_angle` combinations and
accepts an optional non-negative integer `seed`.

```json
{"Rota_Cant": {
  "R_angle": [30, 60],
  "axis": [[0, 1, 0], [0, 0, 1]],
  "C_angle": [45, 90],
  "seed": 12345
}}
```

`Random.num` is a positive integer giving the number of configurations. In
each configuration every nonzero atom independently samples a uniform
direction on the sphere and recovers its original magnitude. Zero moments stay
zero.

```json
{"Random": {"num": 5, "seed": 12345}}
```

`Scale` changes only the magnitude through `m' = (1 + delta)m`. It requires
`0 < pert < 1`, `pert_step > 0`, and an integer `pert / pert_step`. Deltas run
from `-pert` through `+pert`, with negative values first and zero excluded. For
example, `pert=0.25` and `pert_step=0.05` produce ten local variants.

```json
{"Scale": {"pert": 0.25, "pert_step": 0.05}}
```

Each Canting, Rota_Cant, or Random item owns one random-number generator. It is
not reset per atom, variant, INCAR template, or snapshot. Other operation items
do not alter that item's starting moments or random-number consumption.

For example, a Rotation with two angles and two axes, a Canting with two
angles, and the ten-value Scale grid above create `4 + 2 + 10 = 16` independent
perturbed tasks, plus the shared `origin/000000` baseline:

```text
03.spin/.../origin/000000/INCAR
03.spin/.../000-rotation/R1/INCAR
03.spin/.../000-rotation/R4/INCAR
03.spin/.../001-canting/C2/INCAR
03.spin/.../002-scale/S10/INCAR
```

## Writing `machine.json`

`spin_init` accepts the nested `fp` schema used by DP-GEN initial-data
workflows and remains compatible with the legacy flat `fp_*` form:

```json
{
  "api_version": "1.0",
  "fp": {
    "machine": {
      "batch_type": "Slurm",
      "context_type": "local",
      "local_root": "./",
      "remote_root": "/path/to/work"
    },
    "resources": {
      "number_node": 1,
      "cpu_per_node": 32,
      "group_size": 1
    },
    "command": "srun vasp_std",
    "user_forward_files": ["/path/to/KPOINTS"]
  }
}
```

When Stage 2 uses `vasp_std` and Stage 4 uses `vasp_ncl`, add an optional
`spin` entry with the same nested structure. Stage 4 reuses `fp` when `spin` is
absent:

```json
{
  "api_version": "1.0",
  "fp": {
    "machine": {"batch_type": "Slurm", "context_type": "local",
                "local_root": "./", "remote_root": "/path/to/work"},
    "resources": {"number_node": 1, "cpu_per_node": 32, "group_size": 1},
    "command": "srun vasp_std",
    "user_forward_files": ["KPOINTS"]
  },
  "spin": {
    "machine": {"batch_type": "Slurm", "context_type": "local",
                "local_root": "./", "remote_root": "/path/to/work"},
    "resources": {"number_node": 1, "cpu_per_node": 32, "group_size": 1},
    "command": "srun vasp_ncl",
    "user_forward_files": ["KPOINTS"]
  }
}
```

Stage 5 additionally needs a top-level `convert-data` entry:

```json
{
  "convert-data": [{
    "machine": {
      "batch_type": "Slurm",
      "context_type": "local",
      "local_root": "./",
      "remote_root": "/path/to/work"
    },
    "resources": {
      "number_node": 1,
      "cpu_per_node": 2,
      "group_size": 1,
      "queue_name": "partition"
    },
    "command": "nequip-data -m -z 8"
  }]
}
```

Adapt all resources and commands to the cluster. `convert-data` may be one
object or a one-item list. Its terminal simple command must be `nequip-data`.
Environment activation may precede it with `&&`, but no command or pipeline may
follow it, and the entire field must be one line without shell comments. Stage
5 creates remote `out/` and appends `-p data -o out/data.extxyz` when absent.
Explicit values must use those exact paths. DPDispatcher retrieves
`out/data.extxyz`.

Placing `vasp.slurm` in `user_forward_files` only transfers the file into each
task. A command such as `sbatch vasp.slurm` actually invokes sbatch, but ordinary
sbatch returns after enqueueing its child job. DPDispatcher may then retrieve
outputs before VASP has completed. Prefer direct `srun vasp_std` or
`srun vasp_ncl` inside the Slurm job managed by DPDispatcher.

For Bohrium or `DPCloudServerContext`, install the optional dependencies:

```bash
python -m pip install oss2
python -m pip install "dpdispatcher[bohrium]"
```

## Output files

```text
out_dir/
|- param.json
|- 00.scale_pert/
|  `- scale-1.000/
|     |- 000000/POSCAR
|     `- 000001/POSCAR
|- 01.md/
|  |- INCAR
|  |- POTCAR
|  `- scale-1.000/000000/
|     |- POSCAR -> matching structure in 00.scale_pert
|     |- INCAR  -> 01.md/INCAR
|     |- POTCAR -> 01.md/POTCAR
|     |- OUTCAR
|     `- XDATCAR
|- 02.disp/
|  `- scale-1.000/000000/
|     |- 00/POSCAR
|     |- 01/POSCAR
|     `- 02/POSCAR
|- 03.spin/
|  |- POTCAR
|  |- tasks.json
|  `- scale-1.000/000000/00/
|     |- 000-incar/
|     |  |- origin/000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
|     |  |- 000-rotation/R1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
|     |  |- 001-canting/C1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
|     |  `- 002-scale/S1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
|     `- 001-incar/
|        `- ...
`- 04.data/
   |- selection.json
   `- scale-1.000/
      |- data -> 03.spin/scale-1.000/data
      `- out/
         |- data.extxyz
         |- type.raw, type_map.raw, box.raw, coord.raw, energy.raw, force.raw
         |- force_mag.raw, spin.raw, virial.raw
         `- set.000/{box,coord,energy,force,force_mag,spin,virial}.npy
```

Stage 2 always retrieves OUTCAR and XDATCAR. Stage 4 always retrieves OUTCAR
and OSZICAR, and submissions containing relaxation tasks also retrieve CONTCAR.
User backward files are appended. Stage 4 POSCAR/POTCAR entries are real
relative symlinks; each magnetic INCAR is an independent regular file.

The tree above shows list-form `spin_incar`. String input omits the
`###-incar` layer but retains the `origin/000000` baseline and operation-group
layers. Saved manifests from earlier layouts can still be submitted with
`spin_action="run"`; existing trees are not automatically renamed or moved.

See the [Output Directory Reference](../../SPIN_INIT_OUTPUT_STRUCTURE_EN.md)
for complete numbering, symlink, and provenance details.

## Running stages separately

To create Stage 4, inspect its INCARs, and submit later:

```bash
# PARAM: stages=[4], spin_action="make"
dpgen spin_init spin-init.json

# After inspecting 03.spin, change spin_action to "run"
dpgen spin_init spin-init.json machine.json
```

Stage 4 does not read or synchronize `md_incar`, although PARAM retains that
field for the common schema. Existing stage directories are never silently
overwritten. When `make` ran without MACHINE, a later `run` materializes any
missing user-forward symlinks from MACHINE; a conflicting existing file causes
an error instead of being overwritten.

Run-only submission accepts legacy paths saved in `tasks.json`, but does not
regenerate moments, rewrite INCARs, or migrate directories. Use a new `out_dir`
to generate the current layout. A Stage-4-only run still requires matching
`02.disp` snapshots under that output root.

## Stage 5: magnetic result collection and DeepMD export

After every `03.spin` calculation has returned OUTCAR and OSZICAR, set
`"stages": [5]` and run:

```bash
dpgen spin_init spin-init.json machine.json
```

Stage 5 may also be included in `"stages": [1, 2, 3, 4, 5]`. MACHINE must
contain `convert-data`. Stage 5 does not submit VASP; it submits only the data
conversion job. The common PARAM schema still contains the POSCAR, MD INCAR,
and structural perturbation fields, but Stage 5 does not reread those inputs.

Stage 5 first requires every magnetic task in `tasks.json` to have terminated
normally. Every such task enters the exported dataset. `selection.json` records
the source path, scale, and per-scale index for each task; `rejected` is an
empty list.

Within each scale, OUTCAR/OSZICAR pairs are numbered as `OUTCAR-1`,
`OSZICAR-1`, and so on. `convert-data` produces `out/data.extxyz`. `out2npy`
then writes the following files in the same `out/` directory:

```text
type_map.raw
type.raw
box.raw
coord.raw
energy.raw
force.raw
force_mag.raw
spin.raw
virial.raw
set.000/*.npy
```

`energy.npy` has shape `(nframes,)`; other arrays are two-dimensional with one
row per frame. The extxyz and raw intermediates are retained. The number of
frames comes from the returned extxyz and is not assumed to equal the number of
saved tasks. Existing `04.data` is never overwritten. If execution is
interrupted after a remote conversion has been submitted, inspect the remote
job before rerunning to avoid duplicate submission.
