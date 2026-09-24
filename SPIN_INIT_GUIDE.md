# `spin_init` 中文上手手册

最新版入门教程见 [README.md](README.md)。本文保留更详细的配置与排错说明；
字段与完整输入输出说明见
[doc/init/spin-init-usage.md](doc/init/spin-init-usage.md)。
逐层输出目录与文件说明见 [SPIN_INIT_OUTPUT_STRUCTURE.md](SPIN_INIT_OUTPUT_STRUCTURE.md)。

## 1. 当前功能

```text
POSCAR
  → 00.scale_pert：扩胞、缩放、晶胞/原子扰动
  → 01.md：VASP AIMD
  → 02.disp：XDATCAR 的逐帧 POSCAR
  → 03.spin：每个快照的非共线磁性静态或结构/晶格优化 VASP task
  → 04.data：汇总全部已完成任务并导出 DeepMD 磁性数据
```

这是独立命令，不改变 `dpgen init_bulk`。stage 4 已完成输入验证、五种磁矩操作的独立
分组、目录生成、符号链接、dpdispatcher 提交以及 OUTCAR/OSZICAR 回传检查。输入磁矩
不变的基准构型始终放在 `origin/000000/`，与磁扰动模式分组并列。

## 2. 安装与检查

建议使用目标 Linux 集群的 Python 3.9 环境：

```bash
git clone https://github.com/LHCstat/spin_init.git
cd spin_init
conda create -n spin_init python=3.9
conda activate spin_init
python -m pip install .
python -m pip install oss2
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
├── INCAR.spin.1
├── INCAR.spin.2
├── POTCAR
├── KPOINTS
├── spin-init.json
└── machine.json
```

完整流程有两类 INCAR：`INCAR.md` 只供 AIMD 使用；stage 4 可以使用一个或多个
`INCAR.spin.*` 初始磁矩模板。两类文件不会相互覆盖。

每个 stage-4 INCAR 中必须有正确的 `MAGMOM` 和 `M_CONSTR` 标签。每个原子写三个
分量，且两行数值相同。例如二原子体系的静态模板：

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

若 Stage 4 需要结构/晶格优化，把静态参数替换为例如 `NSW=50`、`IBRION=2`，并设置
`ISIF=3`（允许晶格变化）或 `ISIF=2`（固定晶格只优化位置），按体系设置 `EDIFFG`。
程序保留这些参数，不会强制改回静态值。优化任务自动回传 CONTCAR；结果检查将正常
结束与结构收敛分开，未发现结构收敛标志时告警。初始 POSCAR 链接与 `02.disp` 源文件
不会被最终结构覆盖。详见
[Stage 4 结构与晶格优化](doc/init/spin-init-usage.md#stage-4-结构与晶格优化)。

`spin-init.json`：

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
  "md_nstep": 3,
  "spin_incar": ["./INCAR.spin.1", "./INCAR.spin.2"],
  "pert_spin": [
    {"Rotation": {"angle": 45, "axis": [0, 0, 1]}},
    {"Canting": {"angle": [30, 60], "seed": 12345}},
    {"Scale": {"pert": 0.10, "pert_step": 0.10}}
  ],
  "spin_action": "make_run",
  "potcars": ["./POTCAR"]
}
```

`spin_incar` 也可以继续写成单个字符串，此时不增加 `###-incar` 索引层。写成列表时，
列表必须非空、路径不能重复，并按输入顺序增加 `000-incar`、`001-incar`……目录。每个
模板都使用自己的初始 `MAGMOM/M_CONSTR`，然后执行同一套 `pert_spin` 操作。随机操作
共用每个随机操作各自的连续 RNG 序列，遍历顺序为 snapshot → INCAR 输入顺序 →
扰动字段输入顺序 → 局部变体，因此
不同 INCAR 的随机结果不同，但相同 seed 和输入顺序仍可完整复现。

`pert_spin` 是独立操作组列表：每项必须且只能写一个模式；同一模式可以重复。
每个字段从模板的初始磁矩出发，不使用其他字段的结果；只有同一字段内部的参数列表
保留笛卡尔组合。上例产生 `1 + 2 + 2 = 5` 个独立扰动构型，另有共享基准 `origin/000000`。

输出按模式和从零开始的字段输入序号分组：`000-rotation/R1`、`001-canting/C1`、
`001-canting/C2`、`002-scale/S1`、`002-scale/S2`。如果再次输入 Rotation，例如位于
列表第 4 项，则创建独立的 `003-rotation`。分开写 Rotation 和 Canting 不再串联；
需要先旋转再 canting 时使用 `Rota_Cant`。

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
局部变体、原子、后续 INCAR 模板和 snapshot 的稳定顺序连续推进；相同 seed、输入和任务
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
      "queue_name": "partition"
    },
    "command": "nequip-data -m -z 8"
  }]
}
```

旧式 `fp_machine`、`fp_resources`、`fp_command` 等扁平写法也兼容。若 AIMD 使用
`vasp_std`、磁性计算使用 `vasp_ncl`，可在同一 JSON 中再加一个结构完全相同的 `spin`
项，并把它的 command 写为 `srun vasp_ncl`。没有 `spin` 项时，stage 4 自动复用 `fp`。
Stage 5 使用独立的 `convert-data` 项；它从每个 scale 的 `data/` 生成
`out/data.extxyz`，再由 `out2npy` 写入 raw 和 `out/set.000/*.npy`。`command`
最后一条简单命令必须是 `nequip-data`；程序会创建远端 `out/`，并自动补齐缺少的
`-p data -o out/data.extxyz`。若命令自行写出这两个参数，其值必须完全一致。
可在前面使用 `source ... &&` 激活环境，但 `nequip-data` 后不要再接其他命令或管道。
整个 `command` 必须写在一行内，且不能包含 shell 注释。

将 `vasp.slurm` 加入 `user_forward_files` 只会把它传入 task；command 写成
`sbatch vasp.slurm` 则会调用 sbatch 提交该脚本。这不等于 DPDispatcher 会跟踪脚本中的
VASP：普通 sbatch 在子作业入队后就返回，DPDispatcher 可能提前判定命令完成并回传
输出。推荐让 DPDispatcher 管理 Slurm 作业，在其中直接执行 `srun vasp_std` 或
`srun vasp_ncl`，不要把普通异步 sbatch 当作可靠的完成检查方式。

## 5. 分阶段运行

推荐第一次在集群逐阶段检查：

```bash
# stages=[1]
dpgen spin_init spin-init.json machine.json

# stages=[2]：不传 MACHINE 只建目录；传 MACHINE 会提交 AIMD
dpgen spin_init spin-init.json machine.json

# stages=[3]：本地检查 OUTCAR/XDATCAR 并生成 02.disp
dpgen spin_init spin-init.json

# stages=[4], spin_action="make"：只建立磁性任务
dpgen spin_init spin-init.json

# stages=[4], spin_action="run"：提交已有 03.spin
dpgen spin_init spin-init.json machine.json

# stages=[5]：汇总已完成的磁性任务并导出 04.data
dpgen spin_init spin-init.json machine.json
```

也可使用 `stages=[1,2,3,4,5]` 与 `spin_action="make_run"` 一次执行。stage 4 单独运行时
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
├── 03.spin/scale-1.000/000000/00/
│   ├── 000-incar/
│   │   ├── origin/000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│   │   ├── 000-rotation/R1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│   │   ├── 001-canting/C1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│   │   └── 002-scale/S1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│   └── 001-incar/
│       ├── origin/000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       ├── 000-rotation/R1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       ├── 001-canting/C1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       └── 002-scale/S1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
└── 04.data/
    ├── selection.json
    └── scale-1.000/
        ├── data -> 03.spin/scale-1.000/data
        └── out/
            ├── data.extxyz
            ├── {box,coord,energy,force,force_mag,spin,virial}.raw
            ├── type.raw、type_map.raw
            └── set.000/{box,coord,energy,force,force_mag,spin,virial}.npy
```

`01.md` 和 `03.spin` 的 POSCAR/POTCAR 都使用真实相对 symbolic link。stage-4 INCAR
是普通独立文件，每个局部变体目录写入对应的 MAGMOM/M_CONSTR。若 `spin_incar` 是单个
字符串，则不会出现 `###-incar` 层；磁扰动的模式分组层仍按新规则生成。

旧版的 `incar-000` 索引层、`Rotation-000` 等模式分组、直接位于 `000000/` 的基准
任务和未分组扰动任务，仍可使用 `spin_action="run"` 按清单提交；程序不会迁移或重新
扰动已有 `03.spin`。要生成含 `origin/000000/` 的
新布局，请使用新的 `out_dir`；若只运行 Stage 4，
需要先在该输出根目录下准备好对应的 `02.disp` 快照。

stage 2 固定回传 OUTCAR 和 XDATCAR；stage 4 固定回传 OUTCAR 和 OSZICAR，包含优化
任务时这批提交还自动回传 CONTCAR。以上树形结构是静态输出示例。

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

## 8. Stage 5：磁性数据汇总与导出

`03.spin` 中的 VASP 任务全部正常结束且 OUTCAR、OSZICAR 已回传后，
将 PARAM 中的 `stages` 改为 `[5]`，执行：

```bash
dpgen spin_init spin-init.json machine.json
```

也可在完整流程中使用 `stages=[1,2,3,4,5]`。Stage 5 不再次运行 VASP，
但需要 MACHINE 中的 `convert-data` 配置提交 `nequip-data`。Stage 5 不进行 RMSE
筛选；每个 scale 的全部已完成 OUTCAR/OSZICAR 按 `OUTCAR-1`/`OSZICAR-1` 编号收集，
并在 `04.data/scale-*/data` 建立指向 `03.spin/scale-*/data` 的相对链接。
`convert-data` 回传 `04.data/scale-*/out/data.extxyz`；同一个 `out/` 下保留
`*.raw` 和 `set.000/*.npy`。`energy.npy` 为一维；`04.data/selection.json` 记录编号、来源任务和 scale，`rejected` 固定为空。
目前转换提交使用临时工作目录；若在提交后中断，重新执行不会自动恢复原远端任务。
重跑前先确认原任务已结束，避免重复提交。
更多字段与版本约定见 [使用说明](doc/init/spin-init-usage.md#stage-5磁性任务汇总与-deepmd-数据导出)。
