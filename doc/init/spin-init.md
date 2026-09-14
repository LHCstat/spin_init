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

The full workflow uses two separate INCAR files. `md_incar` is the AIMD INCAR
for stage 2. `spin_incar` is the static noncollinear template for stage 4. The
stage-4 template must contain equal `MAGMOM` and `M_CONSTR` vectors, with three
finite Cartesian components for every atom in the POSCAR. It must also enable
`LNONCOLLINEAR` (or `LSORBIT`) and describe a static calculation (`NSW=0` and
`IBRION=-1`).

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
  "pert_spin": [{
    "Canting": {
      "angle": [30, 60],
      "Rcut": [0.4, 0.5],
      "direction": [[1, 0, 0], [0, 1, 0]]
    }
  }],
  "spin_action": "make_run",
  "potcars": ["POTCAR"]
}
```

`spin_action` controls stage 4:

- `make`: create `03.spin` only;
- `run`: submit an existing `03.spin` (a MACHINE file is required);
- `make_run`: create the tasks and submit them when MACHINE is supplied, or
  only create them when MACHINE is omitted.

`pert_spin` currently supports the `Canting` mode. Scalar values or lists are
accepted for `angle` (degrees) and `Rcut`; `direction` accepts one Cartesian
3-vector, a list of 3-vectors, or may be omitted. Parameter lists are expanded
as a Cartesian product and named `C1`, `C2`, and so on. The original moments
remain available as `000000`.

For each nonzero initial moment **a**, Canting projects `direction` into the
plane normal to **a**, rotates that projected direction counterclockwise by
`angle` using the right-hand rule around **a**, and gives it length `Rcut`. The
component along **a** is shortened so that the final moment has the same norm as
the original. Zero moments remain zero. If `direction` is omitted or its
projection vanishes, the Cartesian x/y/z axis least parallel to **a** is chosen
deterministically. `Rcut` may not exceed the magnitude of any affected nonzero
moment.

`spin_pert_numb` is an internal compatibility field and should be left at its
default zero. Rotation, combined Rotation/Canting, Random, and Scale modes are
not implemented yet and are rejected rather than silently ignored.

The output root is `out_dir` exactly; no suffix is appended. Each stage refuses
to overwrite its existing directory. There is no `sys-*` layer because this
workflow accepts one initial POSCAR.

```text
00.scale_pert/scale-1.000/000000/POSCAR
01.md/scale-1.000/000000/{POSCAR,INCAR,POTCAR,OUTCAR,XDATCAR}
02.disp/scale-1.000/000000/00/POSCAR
03.spin/scale-1.000/000000/00/000000/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
03.spin/scale-1.000/000000/00/C1/{POSCAR,INCAR,POTCAR,OUTCAR,OSZICAR}
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
