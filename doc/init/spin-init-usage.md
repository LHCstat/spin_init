# `spin_init` 使用说明

## 作用与流程

`spin_init` 从一个已有 POSCAR 出发，生成结构扰动后的 VASP AIMD 任务，拆分 XDATCAR
轨迹，并为每个轨迹快照建立非共线磁性静态计算：

```text
POSCAR
  → 扩胞、scale、box/atom perturbation
  → VASP AIMD
  → XDATCAR → 独立 POSCAR snapshots
  → 非共线磁性静态 VASP task
```

命令为：

```bash
dpgen spin_init PARAM [MACHINE]
```

四个 stage 分别是：

| stage | 作用 | 输出目录 |
| --- | --- | --- |
| 1 | 扩胞、缩放、结构扰动 | `00.scale_pert` |
| 2 | 建立 AIMD task；提供 MACHINE 时提交 | `01.md` |
| 3 | 检查 OUTCAR、解析 XDATCAR、导出快照 | `02.disp` |
| 4 | 为每个快照建立/提交磁性静态计算 | `03.spin` |

## 输入文件

完整流程需要两个不同的 INCAR：

| 文件 | 作用 |
| --- | --- |
| `POSCAR` | 初始结构 |
| `INCAR.md` | stage 2 的 AIMD 参数 |
| `INCAR.spin` | stage 4 的非共线磁性静态计算模板 |
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

这是独立的 stage-4 模板，不会替代或修改 `INCAR.md`。目前要求：

- 使用正确标签 `MAGMOM`，不是 `MAMGOM`；
- 同时存在 `MAGMOM` 和 `M_CONSTR`，且两者数值完全相同；
- 每个 POSCAR 原子对应三个笛卡尔分量，所以每行总计必须有 `3 × 原子数` 个数；
- 启用 `LNONCOLLINEAR = .TRUE.`（启用 `LSORBIT` 也满足非共线要求）；
- 当前为静态计算，要求 `NSW = 0`、`IBRION = -1`。

二原子示例：

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

### `spin-init.json` 编写方式

完整四阶段示例：

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
| `spin_incar` | 磁性静态 INCAR 模板路径 |
| `pert_spin` | 有序磁矩操作列表；支持 `Rotation`、`Canting`、`Rota_Cant`、`Random`、`Scale` |
| `spin_pert_numb` | 内部兼容字段；用户应省略或保持为 `0` |
| `spin_action` | `make`、`run` 或 `make_run` |

`spin_action=make` 只生成 `03.spin`；`run` 只提交已存在的 `03.spin`，并要求提供
MACHINE；`make_run` 在提供 MACHINE 时生成后提交，未提供 MACHINE 时只生成。

每个快照保留输入磁矩不变的基准构型 `000000`。`pert_spin` 每项必须且只能包含一个
模式，按列表顺序执行，允许重复模式。每一步都对当前全部分支做笛卡尔展开，只输出
最终叶子。局部名称分别为 `R#`、`C#`、`RC#`、`Rand#`、`S#`，组合名称按执行顺序
用 `-` 连接，例如 `R1-C2-S1`。

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

每个 Canting、Rota_Cant 或 Random 操作各自拥有一个 RNG。它不会为每个原子、变体、
父分支或 snapshot 重新设 seed，而是按稳定顺序连续推进。

例如 Rotation 有 2 个 angle 和 2 个 axis、Canting 有 2 个 angle、Scale 使用上述
10 个增量时，共产生 `4 × 2 × 10 = 80` 个最终扰动构型，另加 `000000`。代表路径为：

```text
03.spin/.../000000/INCAR
03.spin/.../R1-C1-S1/INCAR
03.spin/.../R4-C2-S10/INCAR
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

`vasp.slurm` 放入 `user_forward_files` 只表示把文件传到 task 目录。只有 command 明确
写成例如 `sbatch vasp.slurm` 时才会通过该脚本提交；`srun vasp_ncl` 会直接运行程序。

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
└── 03.spin/
    ├── POTCAR
    ├── tasks.json
    └── scale-1.000/000000/00/
        ├── 000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
        ├── R1-C1-S1/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
        └── R1-C2-S2/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
```

stage 2 固定回传 OUTCAR 和 XDATCAR；stage 4 固定回传 OUTCAR 和 OSZICAR。用户配置的
backward files 会在此基础上追加。`03.spin` 中 POSCAR/POTCAR 必须是真实相对符号链接，
INCAR 则是每个磁构型自己的普通文件。

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

## 当前未实现内容

- 磁矩 RMSE 筛选；
- DeepMD 的 `spin.npy`、`spin_force.npy` 数据转换。
