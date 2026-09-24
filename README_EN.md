# `spin_init` User Guide

English | [Chinese](README.md)

`spin_init` starts from an existing POSCAR, generates structurally perturbed
configurations, runs VASP AIMD, exports every XDATCAR frame as an independent
POSCAR, creates noncollinear magnetic calculations for those snapshots, and
finally produces a DeepMD magnetic dataset.

```text
POSCAR -> 00.scale_pert -> 01.md (AIMD) -> 02.disp (POSCAR snapshots)
       -> 03.spin (magnetic VASP tasks) -> 04.data (DeepMD data)
```

Run the workflow with:

```bash
dpgen spin_init PARAM [MACHINE]
```

`PARAM` is the workflow parameter JSON. `MACHINE` is the DPDispatcher machine
configuration used to submit VASP and data-conversion jobs. The stages can be
run separately or in one invocation.

## 1. Installation and runtime requirements

Install this repository in a Python 3.9 environment on the target Linux
cluster:

```bash
git clone https://github.com/LHCstat/spin_init.git
cd spin_init
conda create -n spin_init python=3.9
conda activate spin_init
python -m pip install .
python -m pip install oss2
dpgen spin_init -h
```

Submitting calculations also requires a working VASP installation, available
compute resources, and a valid DPDispatcher configuration. For Bohrium or
`DPCloudServerContext`, install the cloud dependencies with:

```bash
python -m pip install "dpdispatcher[bohrium]"
```

Stage 2 and Stage 4 use real symbolic links for task inputs. Run the production
workflow in a Linux environment in which the account can create symlinks.

## 2. Input files

A typical working directory contains:

```text
work/
|- POSCAR            initial structure
|- INCAR.md          AIMD INCAR
|- INCAR.spin        magnetic INCAR template
|- POTCAR            pseudopotential, or provide separate fragments
|- KPOINTS           forwarded through machine.json when required
|- param_spin.json   spin_init parameters
`- machine.json      execution configuration
```

`INCAR.md` and `INCAR.spin` have different roles. The former is used only for
AIMD; the latter provides the initial magnetic moments and VASP settings for
Stage 4. The atom order in POSCAR, the element order in POTCAR, and the magnetic
moment components in the spin INCAR must agree.

### AIMD INCAR example

```text
SYSTEM = spin_init_md
ENCUT = 500
EDIFF = 1E-5
IBRION = 0
NSW = 2
POTIM = 1.0
TEBEG = 300
TEEND = 300
SMASS = 0
ISMEAR = 1
SIGMA = 0.1
```

### Magnetic INCAR example

The following example describes a two-atom static calculation:

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
```

`MAGMOM` and `M_CONSTR` each contain three Cartesian components per atom, in
POSCAR order. To relax atoms or the cell in Stage 4, set `NSW > 0` and provide
appropriate `IBRION`, `ISIF`, and `EDIFFG` values. Relaxation tasks retrieve
`CONTCAR` in addition to the standard outputs.

## 3. Writing `param_spin.json`

This example runs all five stages:

```json
{
  "stages": [1, 2, 3, 4, 5],
  "from_poscar_path": "./POSCAR",
  "out_dir": "./run_spin",
  "super_cell": [1, 1, 1],
  "scale": [1.0],
  "pert_numb": 1,
  "pert_box": 0.03,
  "pert_atom": 0.01,
  "md_incar": "./INCAR.md",
  "md_nstep": 2,
  "spin_incar": "./INCAR.spin",
  "pert_spin": [
    {"Canting": {"angle": 30, "seed": 12345}}
  ],
  "spin_action": "make_run",
  "potcars": ["./POTCAR"]
}
```

| Field | Meaning |
| --- | --- |
| `stages` | Stages to execute, numbered 1 through 5. Use `[5]` to export completed calculations only. |
| `from_poscar_path` | Path to the initial POSCAR. |
| `out_dir` | Output root used exactly as written, with no automatic suffix. A new empty directory is recommended. |
| `super_cell` | Three positive integers, for example `[2, 2, 2]`. |
| `scale` | Positive lattice scale factors, for example `[0.98, 1.0, 1.02]`. |
| `pert_numb` | Number of randomly perturbed structures per scale; `000000` is additionally retained as the unperturbed scaled structure. |
| `pert_box`, `pert_atom` | Cell and atomic-position perturbation amplitudes. |
| `md_incar`, `md_nstep` | AIMD INCAR and expected step count. If `NSW` differs, the INCAR value is used. |
| `spin_incar` | One magnetic INCAR path, or a nonempty ordered list of paths. |
| `pert_spin` | Independent magnetic perturbation groups. An empty list creates only the `origin` baseline. |
| `spin_action` | Stage 4 action: `make`, `run`, or `make_run`. |
| `potcars` | POTCAR fragment paths in POSCAR element order. |

The workflow never overwrites an existing stage directory. Use a new output
root when creating a workflow. If the stage being created already exists as
`00.scale_pert`, `01.md`, `02.disp`, `03.spin`, or `04.data`, the command stops
with an explicit error. To submit an existing Stage 4 tree, use
`"spin_action": "run"`.

### Multiple initial magnetic templates

To apply the same perturbations to several initial magnetic states, provide an
ordered list:

```json
"spin_incar": ["./INCAR.state_1", "./INCAR.state_2"]
```

Each template receives its own baseline and perturbation tasks. List input adds
an `000-incar`, `001-incar`, ... directory layer. Even a one-item list produces
`000-incar`; a single string omits this layer.

### Magnetic perturbation modes

```json
{
  "pert_spin": [
    {"Rotation": {"angle": [45, 90], "axis": [0, 0, 1]}},
    {"Canting": {"angle": [30, 60], "seed": 12345}},
    {"Rota_Cant": {
      "R_angle": 30,
      "axis": [0, 1, 0],
      "C_angle": 45,
      "seed": 12345
    }},
    {"Random": {"num": 3, "seed": 12345}},
    {"Scale": {"pert": 0.10, "pert_step": 0.05}}
  ]
}
```

- `Rotation` rotates every moment around `axis` by `angle` according to the
  right-hand rule while preserving its magnitude.
- `Canting` places every nonzero moment at the requested angle from its input
  direction. Each atom receives an independently sampled azimuth and retains
  its magnitude.
- `Rota_Cant` applies Rotation first and Canting second within one operation.
- `Random` generates `num` configurations with independently sampled spin
  directions and unchanged input magnitudes.
- `Scale` changes only the magnitude through `m' = (1 + delta)m`. Nonzero
  increments run from `-pert` to `+pert` in steps of `pert_step`.

Every item in `pert_spin` starts independently from the original template
moments. Results from one item are not passed to the next item. Parameter lists
inside one item retain their Cartesian products. `Canting`, `Rota_Cant`, and
`Random` accept a `seed` so that their random sequences can be reproduced.

Task groups use the zero-based input index followed by the lowercase mode, for
example `000-rotation`, `001-canting`, and `002-scale`. Their local variants are
named `R#`, `C#`, `RC#`, `Rand#`, and `S#`. Each snapshot and INCAR template
also has an unchanged `origin/000000` task.

## 4. Writing `machine.json`

The workflow accepts the same nested `fp` configuration used by other DP-GEN
initial-data workflows. Add an optional `spin` entry when Stage 4 needs a
different VASP executable. Replace all paths, partitions, resources, and
commands below with values for the target cluster.

```json
{
  "api_version": "1.0",
  "fp": {
    "machine": {
      "batch_type": "Slurm",
      "context_type": "local",
      "local_root": "./",
      "remote_root": "/path/to/dpdispatcher/work"
    },
    "resources": {
      "number_node": 1,
      "cpu_per_node": 32,
      "group_size": 1,
      "queue_name": "your_partition"
    },
    "command": "srun vasp_std",
    "user_forward_files": ["/path/to/KPOINTS"],
    "user_backward_files": []
  },
  "spin": {
    "machine": {
      "batch_type": "Slurm",
      "context_type": "local",
      "local_root": "./",
      "remote_root": "/path/to/dpdispatcher/work"
    },
    "resources": {
      "number_node": 1,
      "cpu_per_node": 32,
      "group_size": 1,
      "queue_name": "your_partition"
    },
    "command": "srun vasp_ncl",
    "user_forward_files": ["/path/to/KPOINTS"],
    "user_backward_files": []
  },
  "convert-data": [{
    "machine": {
      "batch_type": "Slurm",
      "context_type": "local",
      "local_root": "./",
      "remote_root": "/path/to/dpdispatcher/work"
    },
    "resources": {
      "number_node": 1,
      "cpu_per_node": 2,
      "group_size": 1,
      "queue_name": "your_partition"
    },
    "command": "nequip-data -m -z 8"
  }]
}
```

Stage 4 reuses `fp` when `spin` is absent. The legacy flat `fp_*` schema is
also accepted. `convert-data` may be a one-item list or a single object. Its
terminal simple command must be `nequip-data`; environment setup may precede
it with `&&`, but no pipeline or command may follow it. Stage 5 creates remote
`out/` and adds `-p data -o out/data.extxyz` when those options are omitted.
Explicit path options must use those exact values.

Forwarding `vasp.slurm` only places that file in a task directory. It is
submitted only when `command` explicitly runs `sbatch vasp.slurm`. A plain
`sbatch` returns after enqueueing and may let DPDispatcher retrieve outputs too
early. Prefer direct `srun vasp_std` or `srun vasp_ncl` inside the Slurm job
managed by DPDispatcher.

## 5. Running the workflow by stage

For a first run, execute and inspect each stage separately:

| Stage | PARAM setting and command | Check after completion |
| --- | --- | --- |
| 1: structural perturbation | `stages=[1]`; `dpgen spin_init param_spin.json` | `00.scale_pert/scale-*/000000/POSCAR` and perturbed POSCARs |
| 2: AIMD | `stages=[2]`; `dpgen spin_init param_spin.json machine.json` | `OUTCAR` and `XDATCAR` in every `01.md` task |
| 3: snapshots | `stages=[3]`; `dpgen spin_init param_spin.json` | independent `02.disp/.../00/POSCAR` snapshots |
| 4: create tasks | `stages=[4]`, `spin_action="make"`; omit MACHINE | `03.spin/tasks.json` and task INCARs |
| 4: run tasks | `stages=[4]`, `spin_action="run"`; provide MACHINE | `OUTCAR` and `OSZICAR` in every Stage 4 task |
| 5: export data | `stages=[5]`; `dpgen spin_init param_spin.json machine.json` | `selection.json`, `data.extxyz`, raw files, and `set.000` |

Stage 2 without MACHINE creates the AIMD directories without submitting VASP.
`make_run` creates Stage 4 and submits it when MACHINE is supplied, or only
creates tasks when MACHINE is omitted. Stage 3 checks AIMD completion and
XDATCAR parsing independently. Stage 5 requires the `convert-data` machine
entry and all Stage 4 tasks to have completed normally.

## 6. Output overview

```text
run_spin/
|- param.json
|- 00.scale_pert/scale-1.000/000000/POSCAR
|- 01.md/scale-1.000/000000/{POSCAR,INCAR,POTCAR,OUTCAR,XDATCAR}
|- 02.disp/scale-1.000/000000/00/POSCAR
|- 03.spin/
|  |- POTCAR
|  |- tasks.json
|  `- scale-1.000/000000/00/origin/000000/
|     `- {POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
|- 03.spin/scale-1.000/data/{OUTCAR-1,OSZICAR-1,...}
`- 04.data/
   |- selection.json
   `- scale-1.000/
      |- data -> 03.spin/scale-1.000/data
      `- out/
         |- data.extxyz
         |- *.raw
         `- set.000/*.npy
```

Task-level POSCAR and POTCAR files in `01.md` and `03.spin` are real relative
symlinks. Stage 4 perturbation groups are siblings of `origin`, for example
`000-rotation/R1` and `001-canting/C1`. Multiple `spin_incar` templates add an
`000-incar`, `001-incar`, ... layer between each snapshot and its groups.

Stage 5 includes every normally completed task listed in `tasks.json`, numbers
its OUTCAR/OSZICAR pair within the scale, and records its provenance in
`selection.json`. `convert-data` returns `out/data.extxyz`; `out2npy` then
retains the extxyz and writes raw files plus NumPy arrays under `set.000`.
`energy.npy` is one-dimensional.

For complete directory semantics, numbering rules, and file provenance, see
[Output Directory Reference](SPIN_INIT_OUTPUT_STRUCTURE_EN.md). For every
parameter and Stage 4/5 rule, see the
[Detailed Usage Reference](doc/init/spin-init-usage-en.md).

## 7. Troubleshooting

- `dpgen spin_init` is not recognized: verify that this repository is installed
  in the active environment with `python -m pip show dpgen`.
- `numpy.dtype size changed`: NumPy and h5py are binary-incompatible. Reinstall
  compatible versions in the same environment.
- Bohrium reports `name 'oss2' is not defined`: install `oss2`, or install
  `dpdispatcher[bohrium]` to obtain the Bohrium dependencies.
- A stage directory already exists: choose a new empty `out_dir`, or use an
  operation such as `spin_action="run"` that works with the existing stage.
- Stage 3 cannot find XDATCAR: verify that VASP generated it on the compute node
  and DPDispatcher retrieved it into the corresponding `01.md` task.
- Stage 4 reports a magnetic-vector count mismatch: POSCAR must contain the
  same number of atoms represented by the `3 x natoms` values in both `MAGMOM`
  and `M_CONSTR`.
- `WinError 1314` outside Linux indicates insufficient Windows symlink
  privileges. Production symlink behavior must not be replaced with copies.
