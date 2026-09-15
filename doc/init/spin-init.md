## Spin initialization workflow

`spin_init` is an independent VASP workflow. It preserves `init_bulk` and uses
four stages:

1. create scaled and structurally perturbed POSCARs in `00.scale_pert`;
2. create/run VASP AIMD tasks in `01.md` and retrieve `OUTCAR` and `XDATCAR`;
3. validate the AIMD result, parse `XDATCAR` with pymatgen, and export every
   frame as a POSCAR under `02.disp`;
4. create/run one static noncollinear VASP calculation for every exported
   snapshot under `03.spin`.

```bash
dpgen spin_init PARAM [MACHINE]
```

The full workflow uses two INCAR roles. `md_incar` is the AIMD INCAR for stage
2. `spin_incar` is either one static noncollinear template or an ordered list
of templates for stage 4. Every stage-4 template must contain equal `MAGMOM`
and `M_CONSTR` vectors, with three finite Cartesian components for every atom
in the POSCAR. It must also enable `LNONCOLLINEAR` (or `LSORBIT`) and describe a
static calculation (`NSW=0` and `IBRION=-1`).

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
initial magnetic states with the same pipeline, use:

```json
"spin_incar": ["INCAR.state_1", "INCAR.state_2"]
```

A list must be nonempty and must not repeat a resolved path. List entries add
an `incar-000`, `incar-001`, ... directory between the snapshot and magnetic
configuration. A string keeps the original directory layout. Processing order
is snapshot, then INCAR list order, then perturbation branch. All templates
share the provider's continuous random-number sequence, so their stochastic
results differ but remain reproducible for identical seeds and input order.

`spin_action` controls stage 4:

- `make`: create `03.spin` only;
- `run`: submit an existing `03.spin` (a MACHINE file is required);
- `make_run`: create the tasks and submit them when MACHINE is supplied, or
  only create them when MACHINE is omitted.

`pert_spin` is an ordered operation pipeline. Every list item must contain
exactly one mode, items run from first to last, and repeated modes are allowed.
Each operation expands every current branch over its local parameter variants;
only final leaves are written. The example therefore creates
`1 Rotation * 2 Canting * 2 Scale = 4` perturbed configurations, named
`R1-C1-S1` through `R1-C2-S2`, plus the unchanged `000000` baseline.

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
Each stochastic operation owns one generator that advances through parent
branches, local variants, atoms, and later snapshots in stable order. The same
seed, inputs, and task order reproduce the complete sequence; omitted seeds are
nondeterministic. Zero moments always stay zero. The former Canting `Rcut` and
`direction` fields are rejected.

Local names are `R#`, `C#`, `RC#`, `Rand#`, and `S#`; composed names join
these tokens in execution order with `-`. `spin_pert_numb` remains an internal
compatibility field and should be omitted or left at zero.

The output root is `out_dir` exactly; no suffix is appended. Each stage refuses
to overwrite its existing directory. There is no `sys-*` layer because this
workflow accepts one initial POSCAR.

```text
00.scale_pert/scale-1.000/000000/POSCAR
01.md/scale-1.000/000000/{POSCAR,INCAR,POTCAR,OUTCAR,XDATCAR}
02.disp/scale-1.000/000000/00/POSCAR
03.spin/scale-1.000/000000/00/000000/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
03.spin/scale-1.000/000000/00/R1-C1-S1/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
```

With a `spin_incar` list, the last two paths instead become, for example:

```text
03.spin/scale-1.000/000000/00/incar-000/000000/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
03.spin/scale-1.000/000000/00/incar-001/R1-C1-S1/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
```

Task-level POSCAR and POTCAR files in `01.md` and `03.spin` are real relative
symbolic links. Each stage-4 INCAR is a private regular file so a later magnetic
perturbation can change it independently.

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
adds `OUTCAR` and `OSZICAR`. User backward files are appended without
duplicates. A forwarded `vasp.slurm` merely exists in the task directory; it is
submitted only if `command` explicitly invokes it, such as `sbatch vasp.slurm`.

Stage 3 keeps VASP completion checking independent from XDATCAR parsing. Stage
4 similarly checks static VASP termination independently of future magnetic
quality metrics; RMSE filtering and DeepMD spin-data conversion are not yet
implemented.
