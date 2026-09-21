# `spin_init` 使用说明

## 作用与流程

`spin_init` 从一个已有 POSCAR 出发，生成结构扰动后的 VASP AIMD 任务，拆分 XDATCAR
轨迹，并为每个轨迹快照建立非共线磁性静态或结构/晶格优化计算：

```text
POSCAR
  → 扩胞、scale、box/atom perturbation
  → VASP AIMD
  → XDATCAR → 独立 POSCAR snapshots
  → 非共线磁性静态或优化 VASP task
  → 磁矩 RMSE 筛选 → DeepMD 磁性数据
```

命令为：

```bash
dpgen spin_init PARAM [MACHINE]
```

五个 stage 分别是：

| stage | 作用 | 输出目录 |
| --- | --- | --- |
| 1 | 扩胞、缩放、结构扰动 | `00.scale_pert` |
| 2 | 建立 AIMD task；提供 MACHINE 时提交 | `01.md` |
| 3 | 检查 OUTCAR、解析 XDATCAR、导出快照 | `02.disp` |
| 4 | 为每个快照建立/提交磁性静态或优化计算 | `03.spin` |
| 5 | 检查磁性结果、按 RMSE 筛选并导出 DeepMD 数据 | `04.data` |

## 输入文件

完整流程需要两类 INCAR，其中 stage 4 可以提供一个或多个初始磁矩模板：

| 文件 | 作用 |
| --- | --- |
| `POSCAR` | 初始结构 |
| `INCAR.md` | stage 2 的 AIMD 参数 |
| `INCAR.spin` 或多个 `INCAR.state_*` | stage 4 的非共线磁性静态或优化计算模板 |
| `POTCAR` 或多个赝势片段 | VASP 赝势，顺序必须与 POSCAR 一致 |
| `spin-init.json` | 工作流参数 |
| `machine.json` | 自动提交时使用的 dpdispatcher 配置 |
| `KPOINTS` 等 | 通过 machine 的 `user_forward_files` 传入 |

### `INCAR.md` 编写方式

这是 AIMD INCAR，例如最小测试可写成：

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

若 `NSW` 与 PARAM 的 `md_nstep` 不同，程序沿用原 `init_bulk` 逻辑，以 `NSW`
为准。

### `INCAR.spin` 编写方式

这是独立的 stage-4 模板，不会替代或修改 `INCAR.md`。每个模板都要求：

- 使用正确标签 `MAGMOM`，不是 `MAMGOM`；
- 同时存在 `MAGMOM` 和 `M_CONSTR`，且两者数值完全相同；
- 每个 POSCAR 原子对应三个笛卡尔分量，所以每行总计必须有 `3 × 原子数` 个数；
- 启用 `LNONCOLLINEAR = .TRUE.`（启用 `LSORBIT` 也满足非共线要求）；
- `NSW >= 0`，可自行设置 `NSW / IBRION / ISIF`，程序不会强制改成静态计算。

二原子静态计算示例（不是对所有模板的强制设置）：

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

也支持 VASP 的重复写法和反斜杠续行，例如 `6*0.0`。程序读取 `MAGMOM` 后会验证
其原子数，并为每个 task 写出一份独立 INCAR；`M_CONSTR` 始终与该 task 的 `MAGMOM`
同步。

### Stage 4 结构与晶格优化

需要优化时，在上述 `INCAR.spin` 中替换静态参数并加入例如：

```text
NSW = 50
IBRION = 2
ISIF = 3
EDIFFG = -0.02
```

`NSW` 是最大离子步数，不要求实际运行满 50 步。`IBRION=1/2/3` 配合 `NSW>0`
表示常用结构优化；`ISIF=2` 只优化原子位置，`ISIF=3` 允许位置、晶胞形状和体积变化。
这些仅是写法示例，具体收敛参数需要为体系验证。省略 `IBRION` 时，VASP 在 `NSW>0`
下默认使用 `IBRION=0`，不是结构优化，程序不会自动代填 `IBRION=2`。
参见 [VASP IBRION](https://vasp.at/wiki/index.php/IBRION)、
[ISIF](https://vasp.at/wiki/index.php/ISIF) 和 [NSW](https://vasp.at/wiki/index.php/NSW)。

优化任务自动追加回传 `CONTCAR`。多个模板可同时包含静态与优化模板；由于共享提交
使用统一 backward files，只要有一个优化 task，就为这批所有 task 请求回传 CONTCAR，
但只对优化 task 校验最终结构和收敛标志。用户额外 backward files 继续追加并去重。

结果检查分成两个判断：

- 正常结束：OUTCAR 有一个计时结束标志和力输出、OSZICAR 非空。静态 `NSW=0` 保留
  单个力块检查；离子计算允许多个力块，也允许优化提前结束。
- 优化收敛：优化任务还必须有非空、可解析、原子数及元素顺序正确的 CONTCAR；
  晶格和坐标允许变化。OUTCAR 出现 `reached required accuracy - stopping structural
  energy minimisation` 时记录已确认结构收敛；正常退出但没有该标志时明确告警
  “structural convergence was not confirmed”，不宣称优化收敛。

`check_spin_results` 返回正常结束的任务数，包含发出未确认结构收敛告警的任务；这不代表
电子/磁矩已经收敛。`CONTCAR` 是最后一个离子步结构，即使正常结束也可能未优化收敛。
初始 `POSCAR` 仍为指向 `02.disp` 的真实相对符号链接，程序绝不把 CONTCAR 覆盖到
这个链接或其源文件上。参见 [VASP CONTCAR](https://vasp.at/wiki/index.php/CONTCAR)。

### `spin-init.json` 编写方式

前四阶段示例（需要导出数据时可在 `stages` 末尾添加 `5`）：

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

stage-4 新字段：

| 字段 | 含义 |
| --- | --- |
| `spin_incar` | 一个磁性静态/优化 INCAR 路径，或按顺序排列的非空路径列表 |
| `pert_spin` | 按输入顺序输出的独立磁矩扰动组；支持 `Rotation`、`Canting`、`Rota_Cant`、`Random`、`Scale` |
| `spin_pert_numb` | 内部兼容字段；用户应省略或保持为 `0` |
| `spin_action` | `make`、`run` 或 `make_run` |

`spin_action=make` 只生成 `03.spin`；`run` 只提交已存在的 `03.spin`，并要求提供
MACHINE；`make_run` 在提供 MACHINE 时生成后提交，未提供 MACHINE 时只生成。

单个初始 INCAR 沿用字符串写法，不增加 `incar-###` 索引层：

```json
"spin_incar": "./INCAR.spin"
```

若要对多个初始磁矩状态执行同一套扰动，按需要的处理顺序写成列表：

```json
"spin_incar": [
  "./INCAR.state_1",
  "./INCAR.state_2"
]
```

列表必须非空，每项必须是非空字符串，解析后的文件路径不能重复。列表输入会在
snapshot 与磁构型之间增加 `incar-000`、`incar-001`……层；即使列表中只有一个文件也
会增加 `incar-000`。每个文件提供自己的初始 `MAGMOM/M_CONSTR`，随后应用相同的
`pert_spin`。遍历顺序固定为 snapshot → INCAR 输入顺序 → 扰动字段输入顺序 → 局部变体。
所有 INCAR 共用每个随机操作各自的连续 RNG 序列，不会分别重新 seed，所以随机结果
不同但整体可以复现。

每个快照的每个 INCAR 模板保留一个输入磁矩不变的基准构型，放在 `origin/000000/`。
`origin` 是新增的分类层，与 `Rotation-000` 等模式分组并列；列表形式的多个 INCAR
分别使用 `incar-###/origin/000000/`。上游结构扰动编号和 AIMD 快照编号不变。
`pert_spin` 每项必须且只能包含一个
模式，允许重复模式。不同字段独立作用于模板中的初始磁矩，不会使用前一个字段的
结果；只有同一字段内部的参数列表保留笛卡尔组合。因此不同字段的构型数量相加，
不再相乘。

输出组名为 `<模式>-<从零开始的字段输入序号>`，例如 `Rotation-000`、`Canting-001`、
`Scale-002`。组内局部名称分别为 `R#`、`C#`、`RC#`、`Rand#`、`S#`，例如
`Rotation-000/R1/INCAR`。重复输入同一模式时仍是独立组，不会重名或累积作用。
若需要明确的先旋转再 canting，请使用 `Rota_Cant`，而不是分开写 Rotation 和 Canting。

### Canting 参数和数学定义

`angle` 是扰动后磁矩和初始磁矩之间的夹角，单位为度，取值范围为 `[0, 180]`。
它可以是单个数或列表，每个 angle 生成一个构型。例如 `[30, 60]` 生成 `C1`、`C2`。

对于每个非零初始磁矩 **a**，程序在垂直 **a** 的平面内建立正交基 **e1**、**e2**，
然后为每个原子独立地从 `[0, 2π)` 均匀采样方位角 `phi`：

```text
a' = |a| [cos(angle) a_hat + sin(angle) (cos(phi) e1 + sin(phi) e2)]
```

因此 `a'` 与 `a` 的夹角严格为 `angle`，模长仍为 `|a|`。同一个 INCAR 中不同非零
原子独立采样随机方位角；零磁矩保持不变。
`angle=0` 精确保留原磁矩，`angle=180` 精确反转非零磁矩；这两个端点不消耗随机数。

`seed` 可省略，也可以写成非负整数。相同 seed、输入和任务顺序会复现完整的随机
构型序列；省略 seed 时，每次运行产生新的随机结果。旧的 `Rcut` 和 `direction`
参数不再接受。

### Rotation、Rota_Cant、Random 和 Scale

`Rotation` 的 `angle` 可以是 `[0,360]` 内的一个数或非空列表；`axis` 可以是一个
非零三维向量或向量列表。程序按 `angle × axis` 生成局部变体，并按右手定则使用
Rodrigues 公式绕全局轴旋转所有磁矩。零磁矩不变，非零磁矩模长不变。

```json
{"Rotation": {"angle": [45, 90], "axis": [[0, 0, 1], [0, 1, 0]]}}
```

`Rota_Cant` 是固定顺序的组合操作：先 Rotation，再相对于旋转后的磁矩做 Canting。
它按 `R_angle × axis × C_angle` 的顺序生成所有组合，并支持可选非负整数 `seed`。

```json
{"Rota_Cant": {
  "R_angle": [30, 60],
  "axis": [[0, 1, 0], [0, 0, 1]],
  "C_angle": [45, 90],
  "seed": 12345
}}
```

`Random` 的 `num` 是正整数，表示生成多少组随机方向。每一组中，每个非零原子独立
均匀采样一个球面方向，并乘回该原子的原始磁矩模长；零磁矩保持为零。

```json
{"Random": {"num": 5, "seed": 12345}}
```

`Scale` 只改变模长，使用相对放缩 `m' = (1 + delta) m`。要求 `0 < pert < 1`、
`pert_step > 0`，且 `pert / pert_step` 必须是整数。增量从 `-pert` 到 `+pert`，先负后
正并排除零。例如 `pert=0.25`、`pert_step=0.05` 产生 10 个局部变体。

```json
{"Scale": {"pert": 0.25, "pert_step": 0.05}}
```

每个 Canting、Rota_Cant 或 Random 字段各自拥有一个 RNG。它不会为每个原子、变体、
INCAR 模板或 snapshot 重新设 seed，而是按稳定顺序连续推进。其他字段不会改变该
字段的初始磁矩或随机数消耗次数。

例如 Rotation 有 2 个 angle 和 2 个 axis、Canting 有 2 个 angle、Scale 使用上述
10 个增量时，共产生 `4 + 2 + 10 = 16` 个独立扰动构型，另加共享基准 `origin/000000`。
代表路径为：

```text
03.spin/.../origin/000000/INCAR
03.spin/.../Rotation-000/R1/INCAR
03.spin/.../Rotation-000/R4/INCAR
03.spin/.../Canting-001/C2/INCAR
03.spin/.../Scale-002/S10/INCAR
```

## machine.json 编写方式

`spin_init` 默认沿用 `init_bulk` 的 `fp` 写法，旧式扁平 `fp_*` 写法也继续兼容：

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

如果 stage 2 使用 `vasp_std`，stage 4 使用 `vasp_ncl`，可增加一个与 `fp` 结构完全
相同的可选 `spin` 项。没有 `spin` 项时，stage 4 复用 `fp`：

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

Stage 5 还需要一个独立的 `convert-data` 项，加入 MACHINE 顶层对象：

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

按集群修改资源和命令；`convert-data` 也可写成单个对象。`command` 的最后一条
简单命令必须是 `nequip-data`；前面可用 `source ... &&` 激活环境，但后面不能再接
其他命令或管道。整个字段必须是无 shell 注释的单行命令。程序在远端创建 `out/`，并在命令未指定时自动追加
`-p data -o out/data.extxyz`；若显式指定，两个路径必须完全一致。dpdispatcher
随后回传 `out/data.extxyz`。

`vasp.slurm` 放入 `user_forward_files` 只表示把文件传到 task 目录。command 明确
写成 `sbatch vasp.slurm` 会调用 sbatch 提交该脚本，但普通 sbatch 入队后就返回，
DPDispatcher 不会自动跟踪子作业的 VASP 完成状态，可能提前回传文件并触发结果检查。
推荐由 DPDispatcher 管理 Slurm 作业，在其作业内直接执行 `srun vasp_std` 或
`srun vasp_ncl`，而不是嵌套普通异步 sbatch。

Bohrium/DPCloudServerContext 还需要 dpdispatcher 的 Bohrium 可选依赖：

```bash
python -m pip install "dpdispatcher[bohrium]"
```

## 输出文件

```text
out_dir/
├── param.json
├── 00.scale_pert/
│   └── scale-1.000/
│       ├── 000000/POSCAR
│       └── 000001/POSCAR
├── 01.md/
│   ├── INCAR
│   ├── POTCAR
│   └── scale-1.000/000000/
│       ├── POSCAR -> 00.scale_pert 中对应结构
│       ├── INCAR  -> 01.md/INCAR
│       ├── POTCAR -> 01.md/POTCAR
│       ├── OUTCAR
│       └── XDATCAR
├── 02.disp/
│   └── scale-1.000/000000/
│       ├── 00/POSCAR
│       ├── 01/POSCAR
│       └── 02/POSCAR
├── 03.spin/
│   ├── POTCAR
│   ├── tasks.json
│   └── scale-1.000/000000/00/
│       ├── incar-000/
│       │   ├── origin/000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   ├── Rotation-000/R1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   ├── Canting-001/C1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   └── Scale-002/S1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       └── incar-001/
│           ├── origin/000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│           ├── Rotation-000/R1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│           ├── Canting-001/C1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│           └── Scale-002/S1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
└── 04.data/
    ├── selection.json
    └── scale-1.000/
        ├── data -> 03.spin/scale-1.000/data
        └── out/
            ├── data.extxyz
            ├── type.raw、type_map.raw、box.raw、coord.raw、energy.raw、force.raw
            ├── force_mag.raw、spin.raw、virial.raw
            └── set/{box,coord,energy,force,force_mag,spin,virial}.npy
```

stage 2 固定回传 OUTCAR 和 XDATCAR；stage 4 固定回传 OUTCAR 和 OSZICAR，包含优化
任务的提交还自动回传 CONTCAR。用户配置的
backward files 会在此基础上追加。`03.spin` 中 POSCAR/POTCAR 必须是真实相对符号链接，
INCAR 则是每个磁构型自己的普通文件。上图展示列表写法；字符串写法不含
`incar-###` 层，但同样具有 `origin/000000/` 基准目录和磁扰动模式分组层。

旧版清单中直接位于 `000000/` 的基准任务仍可用 `spin_action="run"` 提交。
程序不会自动搬移旧任务或改写旧清单；新布局需要在新的输出目录中生成。

## 分阶段运行

例如先生成 stage-4 输入、人工检查 INCAR，然后提交：

```bash
# PARAM 中设置 stages=[4], spin_action="make"
dpgen spin_init spin-init.json

# 检查 03.spin 后，把 spin_action 改为 "run"
dpgen spin_init spin-init.json machine.json
```

运行 `stages=[4]` 时不会读取或同步 `md_incar`，但 PARAM 为保持统一格式仍保留该字段。
阶段目录已存在时不会被静默覆盖；已有 `03.spin` 应使用 `spin_action=run`。如果 `make`
时没有 MACHINE，`run` 会按 MACHINE 补建缺少的 user-forward 符号链接；已有同名文件
与配置来源不一致时会报错而不是覆盖。

`run` 仍兼容旧版本保存的未分组磁构型路径（包括串联组合名称），不会迁移或重新
扰动已有 `03.spin`。若要按新规则生成任务，请使用新的 `out_dir`；只运行 Stage 4 时，
该根目录下必须先准备好对应的 `02.disp` 快照。

## Stage 5：磁矩筛选与 DeepMD 数据导出

完成 `03.spin` 全部计算并回传 OUTCAR、OSZICAR 后，将 PARAM 改成
`"stages": [5]`，运行 `dpgen spin_init PARAM MACHINE`；也可以在完整流程
中使用 `"stages": [1, 2, 3, 4, 5]`。MACHINE 必须含 `convert-data`
配置。Stage 5 不重新提交 VASP，但会通过 dpdispatcher 提交转换任务。
参数文件仍须保留统一 schema 规定的 POSCAR、MD INCAR、结构扰动等字段，
但 Stage 5 不重新读取这些输入。

先检查所有磁性任务是否正常结束；任何任务缺失或
未完成都会报出其路径，不会混同为 RMSE 淘汰。然后从各任务 INCAR 读取初始
`MAGMOM`，从最终 OUTCAR 的 x/y/z 磁矩表读取末态磁矩，计算初始非零磁矩原子的
模长 RMSE：`sqrt(mean((|M_initial| - |M_final|)^2))`。超过 `5.0e-3` 的任务
只被排除，不影响其余合格任务，并记入 `selection.json`；没有合格任务则报错。

合格任务的 OUTCAR/OSZICAR 按 scale 编号作为 `convert-data` 的输入，
由该程序生成每个 scale 的 `out/data.extxyz`。随后 `out2npy` 一步在
`out/` 写入 `type_map.raw`、`type.raw`、`box.raw`、`coord.raw`、
`energy.raw`、`force.raw`、`force_mag.raw`、`spin.raw`、`virial.raw`，
并在 `out/set/` 写入相应 `.npy`。`energy.npy` 为一维逐帧数组，
其他 `.npy` 为二维逐帧数组。raw 与 extxyz 均保留。
已有 `04.data` 不会被覆盖，重跑前应先检查并处理该目录。
若 Stage 5 在远端转换任务提交后中断，当前版本不会自动恢复原提交；
再次运行前先检查远端任务状态，以免重复提交。
