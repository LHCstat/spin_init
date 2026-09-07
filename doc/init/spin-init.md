## Spin initialization AIMD snapshots

`spin_init` is an independent VASP workflow that prepares structurally perturbed
AIMD tasks and exports every XDATCAR frame as a standalone POSCAR. It does not
modify `init_bulk` and it does not yet apply spin canting or rotation.

```bash
dpgen spin_init PARAM [MACHINE]
```

The stages are:

1. Create scaled and perturbed structures in `00.scale_pert`.
2. Create AIMD tasks in `01.md`; if `MACHINE` is supplied, submit them with
   dpdispatcher and retrieve both `OUTCAR` and `XDATCAR`.
3. Check VASP completion from `OUTCAR`, parse `XDATCAR` with pymatgen, and write
   snapshots to `02.disp`.

The output root is `out_dir` exactly; no suffix is appended. A stage refuses to
overwrite its own existing output directory, so stages can be run separately.
For example, prepare tasks without submitting them using stages 1 and 2, run the
tasks manually, and later use stage 3 to collect the returned files.

Each scale contains `pert_numb + 1` tasks. Task `000000` is the scaled,
unperturbed structure, while `000001` through `00000N` are perturbed structures.
There is no `sys-*` layer because this workflow accepts one initial POSCAR.

```text
00.scale_pert/scale-1.000/000000/POSCAR
01.md/scale-1.000/000000/{POSCAR,INCAR,POTCAR,OUTCAR,XDATCAR}
02.disp/scale-1.000/000000/00/POSCAR
```

The task-level `POSCAR`, `INCAR`, and `POTCAR` are real relative symbolic links.
The common `01.md/POTCAR` is assembled by concatenating `potcars` in the supplied
order.

An example parameter file is:

```json
{
  "stages": [1, 2, 3],
  "from_poscar_path": "POSCAR",
  "out_dir": ".",
  "super_cell": [1, 1, 1],
  "scale": [1.0],
  "pert_numb": 1,
  "pert_box": 0.03,
  "pert_atom": 0.01,
  "md_incar": "INCAR.md",
  "md_nstep": 3,
  "potcars": ["POTCAR"]
}
```

If `NSW` in `md_incar` differs from `md_nstep`, `spin_init` follows `NSW`, as
`init_bulk` does.

Machine parameters use the normal DP-GEN `fp` configuration. `XDATCAR` is always
added to dpdispatcher's backward files, so it need not be listed in
`user_backward_files`. A forwarded `vasp.slurm` file is only made available in
the task directory. It is submitted only when `fp.command` explicitly requests
that behavior, for example `sbatch vasp.slurm`; a command such as `vasp_std`
runs VASP directly.
