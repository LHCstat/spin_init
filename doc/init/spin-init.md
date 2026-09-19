## Spin initialization workflow

`spin_init` is an independent VASP workflow. It preserves `init_bulk` and uses
five stages:

1. create scaled and structurally perturbed POSCARs in `00.scale_pert`;
2. create/run VASP AIMD tasks in `01.md` and retrieve `OUTCAR` and `XDATCAR`;
3. validate the AIMD result, parse `XDATCAR` with pymatgen, and export every
   frame as a POSCAR under `02.disp`;
4. create/run baseline and independently perturbed noncollinear VASP
   calculations for every exported snapshot under `03.spin`;
5. filter final magnetic results and write DeepMD raw/npy data under
   `04.data/deepmd`.

```bash
dpgen spin_init PARAM [MACHINE]
```

The full workflow uses two INCAR roles. `md_incar` is the AIMD INCAR for stage
2. `spin_incar` is either one noncollinear template or an ordered list
of templates for stage 4. Every stage-4 template must contain equal `MAGMOM`
and `M_CONSTR` vectors, with three finite Cartesian components for every atom
in the POSCAR. It must also enable `LNONCOLLINEAR` (or `LSORBIT`). Both static
calculations and structural/cell relaxations are supported. User `NSW`,
`IBRION`, `ISIF` and other nonmagnetic tags are preserved, not forced to static
values. Standard relaxation modes are `NSW>0` with `IBRION=1`, `2` or `3`;
`ISIF=2` keeps the cell fixed, while `ISIF=3` permits cell shape/volume updates.
See [VASP ISIF documentation](https://vasp.at/wiki/index.php/ISIF).

```json
{
  "stages": [1, 2, 3, 4],
  "from_poscar_path": "POSCAR",
  "out_dir": ".",
  "super_cell": [1, 1, 1],
  "scale": [1.0],
  "pert_numb": 1,
  "pert_box": 0.03,
  "pert_atom": 0.01,
  "md_incar": "INCAR.md",
  "md_nstep": 3,
  "spin_incar": "INCAR.spin",
  "pert_spin": [
    {"Rotation": {"angle": 45, "axis": [0, 0, 1]}},
    {"Canting": {"angle": [30, 60], "seed": 12345}},
    {"Scale": {"pert": 0.10, "pert_step": 0.10}}
  ],
  "spin_action": "make_run",
  "potcars": ["POTCAR"]
}
```

The original string form remains backward compatible. To perturb several
initial magnetic states with the same perturbation groups, use:

```json
"spin_incar": ["INCAR.state_1", "INCAR.state_2"]
```

A list must be nonempty and must not repeat a resolved path. List entries add
an `incar-000`, `incar-001`, ... directory between the snapshot and magnetic
configuration. A string omits this INCAR index layer. Processing order
is snapshot, then INCAR list order, then operation group and local variant. All templates
share each stochastic operation's continuous random-number sequence, so their stochastic
results differ but remain reproducible for identical seeds and input order.

`spin_action` controls stage 4:

- `make`: create `03.spin` only;
- `run`: submit an existing `03.spin` (a MACHINE file is required);
- `make_run`: create the tasks and submit them when MACHINE is supplied, or
  only create them when MACHINE is omitted.

`pert_spin` is a list of independent operation groups. Every item must contain
exactly one mode, and repeated modes are allowed. Each block starts from the
template's original moments, never from another block's result. Cartesian
products apply only to parameter lists inside one block. The example creates
`1 Rotation + 2 Canting + 2 Scale = 5` perturbed configurations, plus one
unchanged `origin/000000` baseline per snapshot/template.

Groups are named by mode and zero-based input index: `Rotation-000`,
`Canting-001`, `Scale-002`. Their local variants are placed in separate task
subdirectories, for example `Rotation-000/R1` or `Canting-001/C2`.
A repeated Rotation at list index 3 gets its own `Rotation-003` group.
Separate Rotation and Canting blocks no longer compose; use `Rota_Cant` when
that explicit combination is required.

Supported modes are:

- `Rotation`: `angle` is one value or list in `[0, 360]`; `axis` is one
  nonzero global 3-vector or a list of vectors. Variants use
  `angle * axis` order and Rodrigues' right-hand-rule rotation.
- `Canting`: `angle` is one value or list in `[0, 180]`. Every nonzero atom
  independently samples an azimuth `phi` in `[0, 2*pi)` and keeps its norm:
  `|m| * (cos(angle) m_hat + sin(angle) (cos(phi) e1 + sin(phi) e2))`.
  Exact 0 and 180 degree endpoints consume no random values.
- `Rota_Cant`: accepts `R_angle`, `axis`, `C_angle`, and optional `seed`.
  Its Cartesian variants use `R_angle * axis * C_angle` order and always
  apply Rotation before Canting.
- `Random`: positive integer `num` generates that many configurations. Every
  nonzero atom receives an independent direction uniformly sampled on the
  sphere while retaining its input norm.
- `Scale`: `0 < pert < 1`, positive `pert_step`, and integer
  `pert / pert_step`. It emits nonzero deltas from `-pert` to `+pert` and
  applies `m' = (1 + delta) m`; negative variants precede positive variants.

Canting, Rota_Cant, and Random accept an optional non-negative integer `seed`.
Each stochastic operation owns one generator that advances through local
variants, atoms, later INCAR templates, and snapshots in stable order. The same
seed, inputs, and task order reproduce the complete sequence; omitted seeds are
nondeterministic. Zero moments always stay zero. The former Canting `Rcut` and
`direction` fields are rejected.

Local names are `R#`, `C#`, `RC#`, `Rand#`, and `S#` inside the corresponding
operation group. `spin_pert_numb` remains an internal
compatibility field and should be omitted or left at zero.

The output root is `out_dir` exactly; no suffix is appended. Each stage refuses
to overwrite its existing directory. There is no `sys-*` layer because this
workflow accepts one initial POSCAR.

```text
00.scale_pert/scale-1.000/000000/POSCAR
01.md/scale-1.000/000000/{POSCAR,INCAR,POTCAR,OUTCAR,XDATCAR}
02.disp/scale-1.000/000000/00/POSCAR
03.spin/scale-1.000/000000/00/origin/000000/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
03.spin/scale-1.000/000000/00/Rotation-000/R1/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
03.spin/scale-1.000/000000/00/Canting-001/C1/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
03.spin/scale-1.000/000000/00/Scale-002/S1/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
04.data/deepmd/<formula>/{type.raw,coord.raw,spin.raw,spin_force.raw,spin_length.raw,set.000/}
```

With a `spin_incar` list, the stage-4 paths instead become, for example:

```text
03.spin/scale-1.000/000000/00/incar-000/origin/000000/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
03.spin/scale-1.000/000000/00/incar-001/Rotation-000/R1/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
```

Task-level POSCAR and POTCAR files in `01.md` and `03.spin` are real relative
symbolic links. Each stage-4 INCAR is a private regular file so a later magnetic
perturbation can change it independently.

Run-only still accepts saved task paths from earlier versions, including
baselines directly under `000000` and ungrouped composed names. Existing
`03.spin` trees are not migrated or re-perturbed; use a new output root to
generate the independent-group layout with the `origin/000000` baseline.
For stage-4-only generation, that root must already contain the required
`02.disp` snapshots.

Machine parameters use the same nested `fp` or legacy flat `fp_*` format as
`init_bulk`. By default both VASP stages use `fp`. If AIMD uses `vasp_std` while
stage 4 must use `vasp_ncl`, an optional `spin` entry may repeat the same nested
shape and override the stage-4 command:

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

Stage 2 always adds `OUTCAR` and `XDATCAR` to backward files. Stage 4 always
adds `OUTCAR` and `OSZICAR`. If any saved task requests standard relaxation,
`CONTCAR` is also added to the shared backward list for all tasks in that
submission (including static tasks in a mixed submission). User backward files are appended without
duplicates. A forwarded `vasp.slurm` merely exists in the task directory; it is
submitted if `command` explicitly invokes sbatch, such as `sbatch vasp.slurm`.
However, ordinary sbatch returns after enqueueing its child job; DPDispatcher
does not automatically track that child's VASP completion and may retrieve
outputs too early. Prefer direct `srun vasp_std` or `srun vasp_ncl` inside the
Slurm job managed by DPDispatcher rather than nesting asynchronous submissions.

Stage 3 keeps VASP completion checking independent from XDATCAR parsing. Stage
4 checks normal termination independently of structural convergence. Static
tasks require one force block; ionic runs may have multiple force blocks and
need not reach `NSW` because optimization can stop early. Relaxation tasks
require a nonempty, parseable `CONTCAR` with the same atom count/species order
as the input. Lattice and coordinates may change. The initial POSCAR symlink
and its snapshot source are never replaced. A normal exit without VASP's
structural convergence marker emits a warning, not a convergence claim;
`check_spin_results` returns the normally terminated task count. This is not
electronic/magnetic convergence or magnetic-quality validation.

Stage 5 can be run alone with `"stages": [5]` after `03.spin` finishes, or
included in `"stages": [1, 2, 3, 4, 5]`. It reads `03.spin/tasks.json` and
requires every saved task to have normally completed before creating output.
Each task contributes only its final ionic frame. DPData reads its standard
structure, energy, and atomic-force labels; the final OUTCAR `magnetization
(x/y/z)` and OSZICAR `MW_int`/`lambda*MW_perp` blocks supply magnetic labels.
VASP 5 uses `2 * lambda*MW_perp`; VASP 6 uses `lambda*MW_perp`. Both are
multiplied by `|MW_int|` without a sign change. Magnetic vectors are rotated
to DPData's Cartesian basis before writing.

For atoms with nonzero initial INCAR moments, the filter computes
`sqrt(mean((|MAGMOM_initial| - |moment_final|)^2))`. Tasks above `5.0e-3` are
listed in `04.data/deepmd/selection.json` and excluded; malformed or
incomplete tasks raise an error. Accepted frames are grouped by composition
and atom order. Each system directory retains ordinary DPData raw/npy files
and aligned `spin.raw`, `spin_force.raw`, `spin_length.raw`, plus matching
`set.000/*.npy` files. `spin` is the unit-vector direction (zero for exactly
zero final moments); `spin_length` stores the original final magnitudes.
`frames.json` maps each dataset row to its source task and OUTCAR frame.
Stage 5 does not submit another VASP calculation or overwrite `03.spin`.
