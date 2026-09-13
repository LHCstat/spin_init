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
  "spin_pert_numb": 0,
  "spin_action": "make_run",
  "potcars": ["./POTCAR"]
}
```

stage-4 新字段：

| 字段 | 含义 |
| --- | --- |
| `spin_incar` | 磁性静态 INCAR 模板路径 |
| `spin_pert_numb` | 额外磁构型数量；扰动规则尚未确定，因此当前必须为 `0` |
| `spin_action` | `make`、`run` 或 `make_run` |

`spin_action=make` 只生成 `03.spin`；`run` 只提交已存在的 `03.spin`，并要求提供
MACHINE；`make_run` 在提供 MACHINE 时生成后提交，未提供 MACHINE 时只生成。

当前不会伪造旋转或倾斜算法。每个快照只产生一个保持输入磁矩不变的基准构型
`000000`。以后实现扰动算法时可在同一级加入 PDF 中建议的 `C1`、`R1` 等名称。

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
    └── scale-1.000/000000/00/000000/
        ├── POSCAR -> 02.disp 中对应快照
        ├── POTCAR -> 03.spin/POTCAR
        ├── INCAR
        ├── OUTCAR
        └── OSZICAR
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

- Canting、Rotation、Rotation & Canting、Random 的物理扰动规则；
- 磁矩 RMSE 筛选；
- DeepMD 的 `spin.npy`、`spin_force.npy` 数据转换。

这些功能需要先确定扰动参数及数学定义，之后可接入已预留的命名与任务展开接口。
