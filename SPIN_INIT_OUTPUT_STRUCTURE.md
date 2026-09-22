# `spin_init` 输出目录结构详解

本文说明 `dpgen spin_init` 生成的全部阶段目录、编号规则、文件来源以及各阶段之间的依赖关系。

以下示例假设参数文件包含：

```json
{
  "out_dir": "./run_spin",
  "scale": [1.0],
  "pert_numb": 1,
  "md_nstep": 3,
  "spin_incar": ["./INCAR.spin.1", "./INCAR.spin.2"],
  "pert_spin": [
    {"Rotation": {"angle": 45, "axis": [0, 0, 1]}},
    {"Canting": {"angle": [30, 60], "seed": 12345}},
    {"Scale": {"pert": 0.10, "pert_step": 0.10}}
  ]
}
```

## 1. 总体数据流

五个阶段按照下面的关系传递结构和计算结果：

```text
输入 POSCAR
    │
    ▼
00.scale_pert                   supercell、scale、晶胞和原子扰动
    │ POSCAR
    ▼
01.md                           每个扰动结构对应一个 VASP AIMD task
    │ OUTCAR + XDATCAR
    ▼
02.disp                         XDATCAR 中的每个 frame 转成独立 POSCAR
    │ POSCAR snapshots
    ▼
03.spin                         对每个 snapshot 的初始磁矩进行独立分组扰动并计算
    │ OUTCAR + OSZICAR
    ▼
04.data                         末态磁矩 RMSE 筛选、extxyz 多帧转换为 raw/npy
```

输出根目录就是 `param.json` 中的 `out_dir`，程序不会在名称后面自动增加后缀。例如：

```json
"out_dir": "./run_spin"
```

对应：

```text
run_spin/
├── param.json
├── 00.scale_pert/
├── 01.md/
├── 02.disp/
├── 03.spin/
└── 04.data/
```

`param.json` 是原始输入参数文件的副本，便于以后追溯计算设置。它不是运行时规范化
参数的完整记录，例如 `md_nstep` 与 INCAR 中 `NSW` 不一致时，实际检查遵循 NSW，
副本不会自动改写该值。分阶段再次运行时，该副本也可能被新输入参数覆盖。

程序不会自动复制 `machine.json` 到输出根目录，建议用户另外保留实际提交用的
machine 配置。

本文目录树中的箭头表示链接的逻辑目标，省略了实际相对链接字符串中的 `../`。
OUTCAR、XDATCAR、OSZICAR 和运行日志只有在执行相应计算并回传后才会出现；
仅运行 `make` 时不会生成这些计算结果文件。

## 2. 目录编号规则

### 2.1 Scale 目录

每个 `scale` 值对应一个目录，固定保留三位小数：

```text
scale-1.000
scale-0.980
scale-1.020
```

例如：

```json
"scale": [0.98, 1.0, 1.02]
```

会在每个阶段产生同名的三个 `scale-*` 分支。

### 2.2 结构扰动编号

结构编号固定为六位数字：

```text
000000    未进行随机 box/atom 扰动的 scaled structure
000001    第 1 个扰动结构
000002    第 2 个扰动结构
...
```

`pert_numb=N` 的含义是生成 `N` 个扰动结构，同时保留一个 `000000` 基准结构。因此，每个 scale 下共有 `N+1` 个结构。

### 2.3 AIMD frame 编号

Stage 3 从零开始为 XDATCAR frame 编号，最少使用两位数字：

```text
00    XDATCAR 中第 1 个 frame
01    XDATCAR 中第 2 个 frame
02    XDATCAR 中第 3 个 frame
...
99
100
```

这里的编号是 AIMD 快照编号，不是结构扰动编号。

### 2.4 磁扰动构型编号

每个 snapshot 的每个初始 INCAR 模板始终包含一个未扰动磁矩基准构型：

```text
origin/000000
```

`origin` 是基准磁矩的分类层，与 `000-rotation` 等分组并列。任务文件在其下的
`000000/`，不在 `origin/` 本身；这个变化不影响 Stage 1–3 的数字编号。

其他构型按 `pert_spin` 字段的模式和输入序号分组：

```text
000-rotation/R1     输入第 1 个字段的第 1 个 Rotation 变体
001-canting/C2      输入第 2 个字段的第 2 个 Canting 变体
002-rota_cant/RC1   输入第 3 个字段的第 1 个 Rota_Cant 变体
003-random/Rand3    输入第 4 个字段的第 3 个 Random 变体
004-scale/S4        输入第 5 个字段的第 4 个 Scale 变体
```

组名后面的序号是从零开始的字段输入索引，最少三位，不是某种模式出现的次数。例如
在第三个字段再次输入 Rotation，会创建 `002-rotation`，不会与 `000-rotation` 重名。

不同字段各自从初始磁矩出发，不会串联或跨字段做笛卡尔积。只有同一个字段内部的
参数列表仍做组合，例如 Rotation 的 `angle × axis`、Rota_Cant 的
`R_angle × axis × C_angle`。`Rota_Cant` 本身仍固定先 Rotation 再 Canting。

每个局部变体目录对应一个独立 VASP task。`000-rotation` 等组目录只是分类层，本身
不直接放置 VASP 输入文件，也不是一个计算 task。

## 3. `00.scale_pert`：结构扰动结果

典型结构如下：

```text
00.scale_pert/
└── scale-1.000/
    ├── 000000/
    │   └── POSCAR
    └── 000001/
        └── POSCAR
```

完整路径可以解释为：

```text
00.scale_pert/scale-1.000/000001/POSCAR
│             │           │      └─ 该结构的独立 POSCAR 文件
│             │           └──────── 第 1 个随机扰动结构
│             └──────────────────── scale = 1.0
└────────────────────────────────── 结构生成阶段
```

这里的 POSCAR 都是普通文件：

- `000000/POSCAR`：完成 supercell 和 scale 后的基准结构；
- `000001...00000N/POSCAR`：在同一 scaled structure 上施加 box perturbation 和 atomic perturbation 后的结构。

Stage 1 完成后，临时的 `POSCAR.supercell`、`POSCAR1.vasp` 等中间文件不会保留在最终目录中。

## 4. `01.md`：VASP AIMD 任务

Stage 2 为 `00.scale_pert` 中的每个结构建立一一对应的 AIMD task：

```text
01.md/
├── INCAR
├── POTCAR
└── scale-1.000/
    ├── 000000/
    │   ├── POSCAR -> 00.scale_pert/scale-1.000/000000/POSCAR
    │   ├── INCAR  -> 01.md/INCAR
    │   ├── POTCAR -> 01.md/POTCAR
    │   ├── OUTCAR
    │   ├── XDATCAR
    │   ├── fp.log
    │   └── [user_forward_files]
    └── 000001/
        ├── POSCAR -> 00.scale_pert/scale-1.000/000001/POSCAR
        ├── INCAR  -> 01.md/INCAR
        ├── POTCAR -> 01.md/POTCAR
        ├── OUTCAR
        ├── XDATCAR
        ├── fp.log
        └── [user_forward_files]
```

其中：

| 文件 | 类型 | 来源或用途 |
|---|---|---|
| `01.md/INCAR` | 普通文件 | `md_incar` 的副本，供所有 AIMD task 共用 |
| `01.md/POTCAR` | 普通文件 | 按 `potcars` 列表顺序拼接得到的公共 POTCAR |
| task `POSCAR` | 相对符号链接 | 指向 `00.scale_pert` 中对应的 POSCAR |
| task `INCAR` | 相对符号链接 | 指向 `01.md/INCAR` |
| task `POTCAR` | 相对符号链接 | 指向 `01.md/POTCAR` |
| `OUTCAR` | 回传文件 | VASP AIMD 主输出，由 DPDispatcher 回传 |
| `XDATCAR` | 回传文件 | VASP AIMD 轨迹，由 DPDispatcher 回传 |
| `fp.log` | 运行日志 | DPDispatcher 配置的标准输出和错误输出日志 |
| 用户 forward 文件 | task 输入 | 来自 `fp.user_forward_files`，例如 `KPOINTS` 或 `vasp.slurm` |

Stage 2 的 POSCAR、INCAR、POTCAR 使用相对链接；额外的 user-forward 文件由
upstream 工具建立到源文件的绝对符号链接。这些额外文件与阶段内的公共 POTCAR
不是同一种链接布局。

Stage 2 固定的 DPDispatcher 文件列表是：

```text
forward_files:
  POSCAR
  INCAR
  POTCAR
  + fp.user_forward_files

backward_files:
  OUTCAR
  XDATCAR
  + fp.user_backward_files
```

VASP 还可能在计算节点产生 CONTCAR、OSZICAR、vasprun.xml、WAVECAR 等其他文件，
但它们不属于 Stage 2 默认回传列表。若确实需要本地保留，应在
`fp.user_backward_files` 中明确追加，而不能仅根据计算节点上的文件推断本地必然存在。

把 `vasp.slurm` 放入 `user_forward_files` 只表示它会出现在 task 目录中。是否执行该脚本由 `machine.json` 的 `command` 决定：

```json
"command": "srun vasp_std"
```

表示直接运行 VASP；而：

```json
"command": "sbatch vasp.slurm"
```

表示调用 sbatch 提交该 Slurm 脚本。但这不等于 DPDispatcher 能跟踪子作业内 VASP
是否完成：普通 sbatch 只等待提交成功，入队后就返回；DPDispatcher 根据命令的返回
状态标记 task 完成，可能在 VASP 尚未运行完时就回传文件并检查输出。

因此上面的 sbatch 写法仅用于解释“实际提交脚本”与“文件仅存在”的区别，不是推荐
运行配置。推荐由 DPDispatcher 管理 Slurm 作业，在该作业内直接执行
`srun vasp_std` 或 `srun vasp_ncl`。如果实验室必须使用额外脚本，需要另外保证命令
等待实际 VASP 完成并正确传递失败状态；当前工作流不提供子作业追踪机制。

## 5. `02.disp`：独立 AIMD 快照

Stage 3 读取每个 AIMD task 的 OUTCAR 和 XDATCAR。只有所有 task 的完整性检查和轨迹解析都通过，程序才会一次性创建 `02.disp`：

```text
02.disp/
└── scale-1.000/
    ├── 000000/
    │   ├── 00/
    │   │   └── POSCAR
    │   ├── 01/
    │   │   └── POSCAR
    │   └── 02/
    │       └── POSCAR
    └── 000001/
        ├── 00/
        │   └── POSCAR
        ├── 01/
        │   └── POSCAR
        └── 02/
            └── POSCAR
```

路径含义如下：

```text
02.disp/scale-1.000/000001/02/POSCAR
│       │           │      │  └─ 独立且合法的 VASP POSCAR
│       │           │      └──── XDATCAR 中第 3 个 frame
│       │           └─────────── 第 1 个结构扰动对应的 AIMD
│       └─────────────────────── scale = 1.0
└─────────────────────────────── 快照收集阶段
```

这些 POSCAR 是由 pymatgen `Xdatcar` 解析轨迹后重新写出的普通文件，不是符号链接。每个文件保留相应 frame 的：

- 晶格；
- 元素及原子顺序；
- 原子坐标。

Stage 3 的两个检查相互独立：

1. 根据 OUTCAR 检查 AIMD 是否完整结束、`TOTAL-FORCE` 块数量是否符合 `md_nstep`；
2. 检查 XDATCAR 是否存在、非空、能够解析，以及所有 frame 的原子数量和元素顺序是否一致。

任何一个 task 失败时都不会留下不完整的 `02.disp`。

## 6. `03.spin`：磁矩扰动与磁性静态/优化任务

### 6.1 使用单个 `spin_incar`

参数写成字符串时：

```json
"spin_incar": "./INCAR.spin"
```

不会增加 INCAR 索引层：

```text
03.spin/
├── POTCAR
├── tasks.json
└── scale-1.000/
    └── 000000/
        └── 00/
            ├── origin/
            │   └── 000000/
            │       ├── POSCAR
            │       ├── POTCAR
            │       ├── INCAR
            │       ├── OUTCAR
            │       ├── OSZICAR
            │       ├── fp.log
            │       └── [user_forward_files]
            ├── 000-rotation/
            │   └── R1/
            │       └── ...
            ├── 001-canting/
            │   ├── C1/
            │   │   └── ...
            │   └── C2/
            │       └── ...
            └── 002-scale/
                ├── S1/
                │   └── ...
                └── S2/
                    └── ...
```

### 6.2 使用多个 `spin_incar`

参数写成列表时：

```json
"spin_incar": [
  "./INCAR.spin.1",
  "./INCAR.spin.2"
]
```

每个输入模板对应一层 `###-incar`：

```text
03.spin/
├── POTCAR
├── tasks.json
└── scale-1.000/
    └── 000000/
        └── 00/
            ├── 000-incar/
            │   ├── origin/
            │   │   └── 000000/
            │   │       ├── POSCAR
            │   │       ├── POTCAR
            │   │       ├── INCAR
            │   │       ├── OUTCAR
            │   │       └── OSZICAR
            │   ├── 000-rotation/
            │   │   └── R1/
            │   │       └── ...
            │   ├── 001-canting/
            │   │   ├── C1/
            │   │   │   └── ...
            │   │   └── C2/
            │   │       └── ...
            │   └── 002-scale/
            │       ├── S1/
            │       │   └── ...
            │       └── S2/
            │           └── ...
            └── 001-incar/
                ├── origin/
                │   └── 000000/
                │       └── ...
                ├── 000-rotation/
                │   └── R1/
                │       └── ...
                ├── 001-canting/
                │   ├── C1/
                │   │   └── ...
                │   └── C2/
                │       └── ...
                └── 002-scale/
                    ├── S1/
                    │   └── ...
                    └── S2/
                        └── ...
```

即使列表中只有一个路径，也会出现 `000-incar`：

```json
"spin_incar": ["./INCAR.spin"]
```

只有字符串写法不会产生该层。

多个模板的遍历顺序固定为：

```text
snapshot → spin_incar 输入顺序 → pert_spin 字段输入顺序 → 局部变体顺序
```

每个模板提供自己的初始 `MAGMOM/M_CONSTR`，每个扰动字段都从这份初始数组出发。
每个随机操作使用各自的连续 RNG 序列，因此不同模板得到不同随机结果；在 seed、
输入和路径顺序相同时，完整结果仍可复现。

例如下面这个路径：

```text
03.spin/scale-1.000/000001/02/001-incar/000-rotation/R2/INCAR
│       │           │      │  │         │            │  └─ 独立磁性 INCAR
│       │           │      │  │         │            └──── 该 Rotation 字段第 2 个变体
│       │           │      │  │         └───────────────── pert_spin 第 1 个字段
│       │           │      │  └─────────────────────────── spin_incar 列表第 2 个模板
│       │           │      └────────────────────────────── 第 3 个 AIMD frame
│       │           └───────────────────────────────────── 第 1 个结构扰动
│       └───────────────────────────────────────────────── scale = 1.0
└───────────────────────────────────────────────────────── 磁性计算阶段
```

`R2` 只有在该 Rotation 字段实际包含至少两个参数变体时才会生成。

### 6.3 Stage 4 文件说明

| 文件 | 类型 | 来源或用途 |
|---|---|---|
| `03.spin/POTCAR` | 普通文件 | 按 `potcars` 顺序拼接的公共 POTCAR |
| `03.spin/tasks.json` | 普通 JSON 文件 | 保存所有可提交 task 的相对路径和稳定顺序 |
| task `POSCAR` | 相对符号链接 | 指向 `02.disp` 中对应 snapshot POSCAR |
| task `POTCAR` | 相对符号链接 | 指向 `03.spin/POTCAR` |
| task `INCAR` | 独立普通文件 | 从相应 `spin_incar` 模板生成，写入该构型的 MAGMOM 和 M_CONSTR |
| `OUTCAR` | 回传文件 | 非共线磁性静态/优化计算的 VASP 主输出 |
| `OSZICAR` | 回传文件 | 磁性计算的迭代/能量摘要 |
| `CONTCAR` | 优化时自动回传 | 最后一个离子步的结构，不一定已优化收敛 |
| `fp.log` | 运行日志 | DPDispatcher 标准输出和错误输出日志 |
| 用户 forward 文件 | 相对符号链接 | Stage 4 machine 配置中的 `user_forward_files` |

Stage 4 固定的 DPDispatcher 文件列表是：

```text
forward_files:
  POSCAR
  INCAR
  POTCAR
  + spin.user_forward_files（存在 spin 配置时）
  或 fp.user_forward_files（复用 fp 时）

backward_files:
  OUTCAR
  OSZICAR
  CONTCAR（包含标准优化任务时自动追加）
  + spin.user_backward_files（存在 spin 配置时）
  或 fp.user_backward_files（复用 fp 时）
```

`tasks.json` 是 `spin_action="run"` 的任务清单。运行阶段从该文件读取准确的 task 路径，
不重新枚举或重新生成磁矩。文件示意如下（这里只列出部分 task）：

```json
{
  "tasks": [
    "scale-1.000/000000/00/000-incar/origin/000000",
    "scale-1.000/000000/00/000-incar/000-rotation/R1",
    "scale-1.000/000000/00/000-incar/001-canting/C1",
    "scale-1.000/000000/00/000-incar/002-scale/S1"
  ]
}
```

清单记录的是最终局部变体目录，不是 `000-rotation` 等分类目录。建议保留原始清单，
不要手工改动路径；这些路径必须对应真实目录并满足规定格式。

新基准 task 的清单路径必须包含 `origin/000000`，不能只写 `origin`。
旧版的 `incar-000` 索引层、`Rotation-000` 等模式分组、不含 `origin` 的基准路径，
以及未分组磁扰动路径仍可按清单运行，不会被自动改名或搬移。新布局只用于新生成的
`03.spin`。

## 7. 完整示例目录树

下面展示 `scale=[1.0]`、`pert_numb=1`、每条 XDATCAR 有三个 frame、两个 `spin_incar`
和三种独立磁扰动字段时的总体结构：

```text
run_spin/
├── param.json
│
├── 00.scale_pert/
│   └── scale-1.000/
│       ├── 000000/POSCAR
│       └── 000001/POSCAR
│
├── 01.md/
│   ├── INCAR
│   ├── POTCAR
│   └── scale-1.000/
│       ├── 000000/
│       │   ├── POSCAR -> 00.scale_pert/.../000000/POSCAR
│       │   ├── INCAR  -> 01.md/INCAR
│       │   ├── POTCAR -> 01.md/POTCAR
│       │   ├── OUTCAR
│       │   └── XDATCAR
│       └── 000001/
│           ├── POSCAR -> 00.scale_pert/.../000001/POSCAR
│           ├── INCAR  -> 01.md/INCAR
│           ├── POTCAR -> 01.md/POTCAR
│           ├── OUTCAR
│           └── XDATCAR
│
├── 02.disp/
│   └── scale-1.000/
│       ├── 000000/
│       │   ├── 00/POSCAR
│       │   ├── 01/POSCAR
│       │   └── 02/POSCAR
│       └── 000001/
│           ├── 00/POSCAR
│           ├── 01/POSCAR
│           └── 02/POSCAR
│
├── 03.spin/
│   ├── POTCAR
│   ├── tasks.json
│   └── scale-1.000/
│       ├── 000000/
│       │   ├── 00/
│       │   │   ├── 000-incar/
│       │   │   │   ├── origin/000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   │   │   ├── 000-rotation/R1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   │   │   ├── 001-canting/
│       │   │   │   │   ├── C1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   │   │   │   └── C2/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   │   │   └── 002-scale/
│       │   │   │       ├── S1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   │   │       └── S2/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
│       │   │   └── 001-incar/
│       │   │       └── ...
│       │   ├── 01/
│       │   │   └── ...
│       │   └── 02/
│       │       └── ...
│       └── 000001/
│           ├── 00/
│           │   └── ...
│           ├── 01/
│           │   └── ...
│           └── 02/
│               └── ...

└── 04.data/
    ├── selection.json
    └── scale-1.000/
        ├── data -> 03.spin/scale-1.000/data
        └── out/
            ├── data.extxyz
            ├── type.raw、type_map.raw
            ├── box.raw、coord.raw、energy.raw、force.raw
            ├── force_mag.raw、spin.raw、virial.raw
            └── set/{box,coord,energy,force,force_mag,spin,virial}.npy
```

## 8. 数量关系

设：

- scale 数量为 `Nscale`；
- `pert_numb = Npert`；
- 每个 AIMD task 含 `Nframe` 个 XDATCAR frame；
- `spin_incar` 模板数量为 `Nincar`；
- `pert_spin` 各独立字段的局部变体数量之和为 `Nspin`，不包括基准构型。

则：

```text
00.scale_pert 中 POSCAR 数量
    = Nscale × (Npert + 1)

01.md 中 AIMD task 数量
    = Nscale × (Npert + 1)

02.disp 中 POSCAR snapshot 数量
    = Nscale × (Npert + 1) × Nframe

03.spin 中磁性 VASP task 数量
    = Nscale × (Npert + 1) × Nframe × Nincar × (Nspin + 1)
```

最后的 `+1` 表示每个 snapshot 和每个初始 INCAR 都会保留 `origin/000000` 基准磁矩构型。

本文参数示例中，`Nspin = 1 + 2 + 2 = 5`，而不是 `1 × 2 × 2`。在三个 frame、
两个 INCAR 模板和两个结构 task 的假设下，共有：

```text
1 × 2 × 3 × 2 × (5 + 1) = 72 个磁性 VASP task
```

如果不同 AIMD task 的 frame 数量不同，则 `02.disp` 和 `03.spin` 的实际数量应分别对每个 task 的 frame 数量求和。
`04.data` 的实际帧数由 `convert-data` 回传的 `data.extxyz` 决定；RMSE
按每个 task 的末态磁矩筛选整个 task，不会假定每个合格 task 恰好只贡献一帧。

## 9. 分阶段运行时的目录依赖

| 要运行的 Stage | 必须已有的输入或目录 | 新生成或使用的目录 |
|---|---|---|
| Stage 1 | 初始 POSCAR | 新建 `00.scale_pert` |
| Stage 2 | 完整 `00.scale_pert`、`md_incar`、POTCAR | 新建 `01.md`；传 MACHINE 时提交任务 |
| Stage 3 | `01.md` 中所有 task 的 OUTCAR 和 XDATCAR | 新建 `02.disp` |
| Stage 4 `make` | 完整 `02.disp`、一个或多个 `spin_incar`、POTCAR | 新建 `03.spin` |
| Stage 4 `run` | 已存在且包含有效 `tasks.json` 的 `03.spin` | 提交已有 task，不重建目录 |
| Stage 4 `make_run` | 与 `make` 相同 | 新建后提交；未传 MACHINE 时只建立目录 |
| Stage 5 | `03.spin/tasks.json`、全部 task 正常结束且 OUTCAR/OSZICAR 已回传、MACHINE 含 `convert-data` | 新建 `04.data/scale-*/out`；提交转换任务，不提交 VASP |

各阶段不会静默覆盖同名阶段目录。例如 `02.disp` 已存在时再次执行 Stage 3 会明确报错；`03.spin` 已存在时应使用 `spin_action="run"` 提交已有任务。

`run` 仍支持旧版本 `tasks.json` 中的 `incar-000` 索引层、`Rotation-000` 等模式分组，
以及未分组磁构型路径（包括旧的 `R1-C1-S1` 等名称）。
它不会改变已有 INCAR、迁移目录或按新参数重新扰动。若要生成独立字段分组的新布局，
请指定新的 `out_dir`，而不是在旧 `03.spin` 上重复 `make`。如果只运行 Stage 4，
新输出根目录下也必须先准备好对应的 `02.disp` 快照。

## 10. 检查输出是否正常

### 10.1 检查 Stage 1 数量

```bash
find run_spin/00.scale_pert -name POSCAR -type f
```

每个 scale 应有 `pert_numb + 1` 个 POSCAR。

### 10.2 检查 AIMD 回传文件

```bash
find run_spin/01.md -name OUTCAR -type f
find run_spin/01.md -name XDATCAR -type f
```

OUTCAR 和 XDATCAR 的数量都应等于 AIMD task 数量。

### 10.3 检查真实符号链接

Linux 上可以执行：

```bash
test -L run_spin/01.md/scale-1.000/000000/POSCAR
readlink run_spin/01.md/scale-1.000/000000/POSCAR
```

Stage 4 同样可以检查：

```bash
find run_spin/03.spin -name POSCAR -type l | head
find run_spin/03.spin -name POTCAR -type l | head
```

POSCAR/POTCAR 必须是实际 symbolic link，不应被替换成普通副本。

### 10.4 检查 snapshot 数量

```bash
find run_spin/02.disp -name POSCAR -type f | wc -l
```

### 10.5 检查磁性任务清单

```bash
python -m json.tool run_spin/03.spin/tasks.json
```

清单中的每个相对路径都应对应一个真实 task 目录，目录内至少应有：

```text
POSCAR
POTCAR
INCAR
```

计算完成后还应有：

```text
OUTCAR
OSZICAR
```

### 10.6 检查磁性数据导出

```bash
python -m json.tool run_spin/04.data/selection.json
find run_spin/04.data -path '*/out/data.extxyz' -type f
find run_spin/04.data -path '*/out/set/spin.npy' -type f
```

`selection.json` 的 `selected` 为合格 task，`rejected` 为 RMSE 超阈值 task；
每个 scale 下的 `out/set/*.npy` 应具有相同的第一维；`energy.npy` 为一维。

## 11. 重要说明

- `00.scale_pert`、`01.md`、`02.disp`、`03.spin`、`04.data` 的数字前缀表示数据处理顺序，不是计算编号。
- Stage 2 的 `000000` 表示未随机扰动的结构，不代表第零个 AIMD frame。
- Stage 3 的 `00` 表示 XDATCAR 的第一个 frame。
- Stage 4 的 `origin/000000` 表示未进行磁矩扰动的基准磁构型。
- 单个字符串形式的 `spin_incar` 不增加 INCAR 索引层；列表形式始终增加 `###-incar` 层。
- 不同 `pert_spin` 字段独立作用于初始磁矩；`<字段序号>-<小写模式>` 是分类层，局部变体才是计算 task。
- 同一字段内部的参数组合仍保留；只有 `Rota_Cant` 明确组合 Rotation 与 Canting。
- `01.md` 固定回传 OUTCAR 和 XDATCAR；`03.spin` 固定回传 OUTCAR 和 OSZICAR。
  包含标准优化任务（`NSW>0, IBRION=1/2/3`）时，整批提交还自动追加回传 CONTCAR。
- `03.spin/tasks.json` 是提交已有磁性任务的重要清单，应与 task 目录一起保存。
- 保存或搬移输出时应保留整个阶段关系和符号链接。单独搬移 `01.md` 或 `03.spin`
  可能使 POSCAR 链接失效；指向输出根目录外的 user-forward 文件还需要保留源文件。
- Stage 4 的独立 INCAR 中 `MAGMOM` 和 `M_CONSTR` 始终写成相同的三分量磁矩数组。
- Stage 5 读取 `03.spin/tasks.json` 与已完成的 VASP 输出，通过 MACHINE 的 `convert-data` 任务回传 extxyz，再生成 `04.data/scale-*/out`；不会重新计算 VASP。

## 12. Stage 4 优化时的输出差异

上面的 task 树是静态输出示例。Stage 4 的 INCAR 可以保留用户设置的 `NSW / IBRION /
ISIF`，例如 `NSW=50, IBRION=2, ISIF=3` 允许结构和晶格优化，不增加新的阶段目录。
每个优化 task 在原有文件之外多一个自动回传的普通文件 `CONTCAR`：

```text
03.spin/.../000-rotation/R1/
├── POSCAR -> 对应 02.disp 快照（初始结构，链接保持不变）
├── POTCAR -> 03.spin/POTCAR
├── INCAR
├── OUTCAR
├── OSZICAR
└── CONTCAR   最后一个离子步的结构（可能未收敛）
```

静态与优化模板混合时，共享 backward files 会为所有 task 请求 CONTCAR；只有优化
task 被要求通过最终结构校验。CONTCAR 必须非空、可被 pymatgen Poscar 解析，并与
初始 POSCAR 的原子数/元素顺序一致，晶格和坐标允许变化。它不会替换初始 POSCAR
符号链接，更不会覆盖 `02.disp` 快照。

正常结束检查允许优化的多个离子步和提前结束；结构收敛单独查看 OUTCAR 的 VASP
收敛标志。正常退出但未出现该标志时程序明确告警，而不是宣称结构已收敛。
返回的正常结束 task 数也包含这些告警任务，不代表电子或磁矩收敛。
具体输入及检查规则见
[Stage 4 结构与晶格优化](doc/init/spin-init-usage.md#stage-4-结构与晶格优化)。

## 13. `04.data`：RMSE 筛选与 extxyz 转换

Stage 5 可单独以 `"stages": [5]` 运行，但需要 MACHINE 中的 `convert-data`
配置。它要求 `03.spin/tasks.json` 中的任务正常结束，再对初始非零磁矩原子比较
INCAR/OUTCAR 磁矩模长。RMSE 为 `sqrt(mean((|m_initial|-|m_final|)^2))`，
阈值 `5.0e-3`。不合格任务只记录在 `selection.json`，不会送入转换。

```text
03.spin/scale-1.000/data/
├── OUTCAR-1 -> 合格磁性任务的 OUTCAR
└── OSZICAR-1 -> 同一任务的 OSZICAR
04.data/
├── selection.json
└── scale-1.000/
    ├── data -> 03.spin/scale-1.000/data
    └── out/
        ├── data.extxyz
        ├── type_map.raw、type.raw
        ├── box.raw、coord.raw、energy.raw、force.raw
        ├── force_mag.raw、spin.raw、virial.raw
        └── set/
            ├── box.npy、coord.npy、energy.npy、force.npy
            └── force_mag.npy、spin.npy、virial.npy
```

`convert-data` 在各 scale 的 `data/` 输入上运行。`nequip-data` 必须是 command 的
最后一条简单命令，command 必须为不含 shell 注释的单行命令；程序自动
创建远端 `out/` 并补齐缺少的 `-p data -o out/data.extxyz`，然后回传该文件；
`out2npy` 一步产生同目录的 raw 与 `set/*.npy`。每个 raw 数值文件每帧一行；
`energy.npy` 是形状为 `(帧数,)` 的一维数组。`spin` 由
`spin_length × initial_magmoms` 得到，`force_mag` 对应
`spin_forces_vert`，`virial = -体积 × stress`。

`selection.json` 含 `rmse_limit`、`selected`、`rejected`；合格记录包括
来源 task、scale、该 scale 内的编号和 RMSE。任务缺失或未正常完成会报错，
不会当成 RMSE 淘汰；全部任务被筛除时也不会生成空数据集。
已有 `04.data` 不会被覆盖。
