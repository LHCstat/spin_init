# `spin_init` 中文上手手册

本文给出安装、配置、分阶段运行和排错步骤。字段与完整输入输出说明见
[doc/init/spin-init-usage.md](doc/init/spin-init-usage.md)。

## 1. 当前功能

```text
POSCAR
  → 00.scale_pert：扩胞、缩放、晶胞/原子扰动
  → 01.md：VASP AIMD
  → 02.disp：XDATCAR 的逐帧 POSCAR
  → 03.spin：每个快照的非共线磁性静态 VASP task
```

这是独立命令，不改变 `dpgen init_bulk`。stage 4 已完成输入验证、五种磁矩操作的有序
组合、目录生成、符号链接、dpdispatcher 提交以及 OUTCAR/OSZICAR 回传检查。输入磁矩
不变的基准构型始终命名为 `000000`。

## 2. 安装与检查

建议使用目标 Linux 集群的 Python 3.9 环境：

```bash
git clone https://github.com/LHCstat/spin_init.git
cd spin_init
python3.9 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -e ".[test]"
dpgen spin_init -h
```

Bohrium/DPCloudServerContext 需要额外安装：

```bash
python -m pip install "dpdispatcher[bohrium]"
```

若出现 `numpy.dtype size changed`，说明当前 numpy 与 h5py 二进制不兼容，应在同一环境
重新安装相容版本，而不是修改 `spin_init` 代码。

## 3. 准备输入

```text
work/
├── POSCAR
├── INCAR.md
├── INCAR.spin
├── POTCAR
├── KPOINTS
├── spin-init.json
└── machine.json
```

完整流程有两个 INCAR：`INCAR.md` 只供 AIMD 使用；`INCAR.spin` 只供 stage 4 使用。
两者不会相互覆盖。

`INCAR.spin` 中必须有正确的 `MAGMOM` 和 `M_CONSTR` 标签。每个原子写三个分量，且
两行数值相同。例如二原子体系：

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

`spin-init.json`：

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

`pert_spin` 是有序操作列表：每项必须且只能写一个模式，按列表顺序依次执行；同一模式
可以重复。每一步都会对已有全部分支做笛卡尔展开，只把最终叶子写入 `03.spin`。上例
产生 `1 × 2 × 2 = 4` 个最终构型：`R1-C1-S1` 到 `R1-C2-S2`，另有基准 `000000`。

五种模式如下：

- `Rotation`：`angle` 为 `[0,360]` 内的数或列表，`axis` 为非零三维轴或轴列表；按右手
  定则用 Rodrigues 公式旋转所有磁矩，局部命名为 `R1`、`R2`……。
- `Canting`：`angle` 为 `[0,180]` 内的数或列表；每个非零原子独立随机方位角并保持
  模长，局部命名为 `C1`、`C2`……。0° 和 180° 不消耗随机数。
- `Rota_Cant`：参数为 `R_angle`、`axis`、`C_angle` 和可选 `seed`；固定先 Rotation
  再 Canting，组合顺序是 `R_angle × axis × C_angle`，局部名为 `RC1`……。
- `Random`：正整数 `num` 表示输出组数；每组为每个非零原子独立均匀采样球面方向，
  保持各自模长，局部名为 `Rand1`……。
- `Scale`：要求 `0 < pert < 1`、`pert_step > 0` 且二者之比为整数；按从负到正、排除
  零的相对增量生成 `S1`……，计算式为 `m' = (1 + delta) m`。

Canting、Rota_Cant、Random 的 `seed` 是可选非负整数。每个随机操作拥有一个 RNG，按
父分支、局部变体、原子和后续 snapshot 的稳定顺序连续推进；相同 seed、输入和任务
顺序可复现完整序列。所有模式都让零磁矩保持为零。旧 Canting 参数 `Rcut`、
`direction` 不再接受；`spin_pert_numb` 是内部字段，用户应省略。

## 4. machine.json

普通情况完全沿用 `init_bulk` 的 `fp` 写法：

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
      "queue_name": "partition"
    },
    "command": "srun vasp_std",
    "user_forward_files": ["/path/to/KPOINTS"],
    "user_backward_files": []
  }
}
```

旧式 `fp_machine`、`fp_resources`、`fp_command` 等扁平写法也兼容。若 AIMD 使用
`vasp_std`、磁性计算使用 `vasp_ncl`，可在同一 JSON 中再加一个结构完全相同的 `spin`
项，并把它的 command 写为 `srun vasp_ncl`。没有 `spin` 项时，stage 4 自动复用 `fp`。

将 `vasp.slurm` 加入 `user_forward_files` 只会把它传入 task。只有 command 写成例如
`sbatch vasp.slurm` 时才会执行该脚本；`srun vasp_std` 或 `srun vasp_ncl` 是直接执行
VASP。

## 5. 分阶段运行

推荐第一次在集群逐阶段检查：

```bash
# stages=[1]
dpgen spin_init spin-init.json

# stages=[2]：不传 MACHINE 只建目录；传 MACHINE 会提交 AIMD
dpgen spin_init spin-init.json machine.json

# stages=[3]：本地检查 OUTCAR/XDATCAR 并生成 02.disp
dpgen spin_init spin-init.json

# stages=[4], spin_action="make"：只建立磁性任务
dpgen spin_init spin-init.json

# stages=[4], spin_action="run"：提交已有 03.spin
dpgen spin_init spin-init.json machine.json
```

也可使用 `stages=[1,2,3,4]` 与 `spin_action="make_run"` 一次执行。stage 4 单独运行时
不会读取 `md_incar`，但参数文件仍保留该字段以维持统一 schema。若 `make` 时没有提供
MACHINE，`run` 会根据 MACHINE 安全补建缺少的 KPOINTS 等 user-forward 符号链接；若
task 中已有同名但内容不同的文件则明确报错，不会覆盖。

## 6. 输出目录

```text
run_spin/
├── param.json
├── 00.scale_pert/scale-1.000/000000/POSCAR
├── 01.md/scale-1.000/000000/
│   ├── POSCAR -> 00.scale_pert 中对应结构
│   ├── INCAR  -> 01.md/INCAR
│   ├── POTCAR -> 01.md/POTCAR
│   ├── OUTCAR
│   └── XDATCAR
├── 02.disp/scale-1.000/000000/
│   ├── 00/POSCAR
│   ├── 01/POSCAR
│   └── 02/POSCAR
└── 03.spin/scale-1.000/000000/00/
    ├── 000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
    ├── R1-C1-S1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
    └── R1-C2-S2/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
```

`01.md` 和 `03.spin` 的 POSCAR/POTCAR 都使用真实相对 symbolic link。stage-4 INCAR
是普通独立文件，每个组合目录写入对应的 MAGMOM/M_CONSTR。

stage 2 固定回传 OUTCAR 和 XDATCAR；stage 4 固定回传 OUTCAR 和 OSZICAR。

## 7. 检查与排错

- stage 3 将“OUTCAR 是否完整完成”与“XDATCAR 能否解析”分开检查；所有 task 通过后
  才创建 `02.disp`。
- stage 4 先验证每个 POSCAR 的原子数、元素顺序、MAGMOM/M_CONSTR 长度和相等性，
  再原子化创建整个 `03.spin`。
- `03.spin` 已存在时，`make`/`make_run` 会拒绝覆盖；若要提交已有任务，使用 `run`。
- Windows 的 `WinError 1314` 是本机符号链接权限问题。正式 symlink 测试和计算应在
  Linux 集群执行，不应把链接改成 copy。
- DPCloudServerContext 报 `name 'oss2' is not defined` 时，安装
  `dpdispatcher[bohrium]`。当前代码会在提交前给出明确提示。

## 8. 后续工作

- 实现磁矩 RMSE 筛选；
- 转换 DeepMD 的 `spin.npy`、`spin_force.npy`。
