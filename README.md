# `spin_init` 使用说明

`spin_init` 从已有 POSCAR 出发，生成结构扰动、运行 VASP AIMD、将 XDATCAR 逐帧转为 POSCAR，再为每帧建立非共线磁性计算，最后汇总全部已完成任务并导出 DeepMD 磁性数据。

```text
POSCAR → 00.scale_pert → 01.md (AIMD) → 02.disp (POSCAR snapshots)
       → 03.spin (磁性 VASP 任务) → 04.data (DeepMD 数据)
```

命令：

```bash
dpgen spin_init PARAM [MACHINE]
```

`PARAM` 是工作流参数 JSON；`MACHINE` 是通过 DPDispatcher 提交 VASP 时使用的机器配置 JSON。可以逐阶段运行，也可以在已配置好 VASP 和机器后连续运行。`spin_init` 不改变原有 `dpgen init_bulk`。

## 1. 安装与运行条件

推荐在目标 Linux 集群的 Python 3.9 环境安装本仓库版本，而不是安装 PyPI 上的上游 `dpgen`：

```bash
git clone https://github.com/LHCstat/spin_init.git
cd spin_init
conda create -n spin_init python=3.9
conda activate spin_init
python -m pip install .
dpgen spin_init -h
```

提交计算还需要可用的 VASP、计算资源和 DPDispatcher 配置。使用 Bohrium/DPCloudServerContext 时，另安装 `python -m pip install "dpdispatcher[bohrium]"`。Stage 2 和 4 的任务目录使用真实符号链接；正式运行请使用允许创建 symlink 的 Linux 环境。

## 2. 准备输入文件

一个最小的工作目录可以是：

```text
work/
├── POSCAR            初始结构
├── INCAR.md          AIMD 的 INCAR
├── INCAR.spin        磁性计算的 INCAR 模板
├── POTCAR            赝势；也可用多个片段
├── KPOINTS           如 VASP 命令需要，通过 machine.json 转发
├── param_spin.json   spin_init 参数
└── machine.json      提交 VASP 时使用
```

`INCAR.md` 和 `INCAR.spin` 是两个不同的输入。前者只用于 AIMD；后者提供初始磁矩并用于 Stage 4。`POSCAR` 中的原子顺序、`POTCAR` 的元素顺序以及磁矩分量顺序必须对应。以下 INCAR 仅用于说明字段，VASP 参数应按具体体系调整。

`INCAR.md` 的短 AIMD 示例：

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

`INCAR.spin` 的二原子、静态计算示例：

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

正确的 VASP 标签是 `MAGMOM`，不是 `MAMGOM`。`MAGMOM` 和 `M_CONSTR` 均需为每个原子提供三个笛卡尔分量，初始值相同。程序会为每个磁扰动任务生成独立的 INCAR，并保持这两组磁矩同步。零初始磁矩在扰动后仍为零。若 Stage 4 需要结构或晶格优化，可自行设置 `NSW>0`、`IBRION`、`ISIF`、`EDIFFG`；程序不会强行改为静态计算。优化时会回传 `CONTCAR`，但正常退出不等于确认优化收敛。

## 3. 编写 `param_spin.json`

下面是从结构到数据导出的完整配置示例；它假设 POSCAR 与 INCAR 中的磁矩数量一致：

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

| 字段 | 写法与作用 |
| --- | --- |
| `stages` | 要执行的阶段，编号 1–5；可仅写 `[5]` 导出已完成的计算 |
| `from_poscar_path` | 已有 POSCAR 的路径 |
| `out_dir` | 输出根目录，名称按原值使用，不追加后缀 |
| `super_cell` | 三个正整数，例如 `[2, 2, 2]` |
| `scale` | 正的晶胞缩放因子列表，例如 `[0.98, 1.0, 1.02]` |
| `pert_numb` | 每个 scale 下的随机扰动结构数；另保留一个未随机扰动的 `000000` |
| `pert_box`、`pert_atom` | 晶胞与原子位置扰动幅度；结构扰动沿用 `init_bulk` 算法 |
| `md_incar`、`md_nstep` | AIMD INCAR 和预期步数；若 `NSW` 不同，以 INCAR 的 `NSW` 为准 |
| `spin_incar` | 一个磁性 INCAR 路径，或多个路径组成的非空列表 |
| `pert_spin` | 磁矩扰动组列表，每项只能指定一个模式；也可为空列表，只生成 `origin` 基准任务 |
| `spin_action` | Stage 4 的 `make`、`run` 或 `make_run` |
| `potcars` | 按 POSCAR 元素顺序排列的 POTCAR 片段路径列表 |

有多个初始磁矩方案时，将 `spin_incar` 写为 `"spin_incar": ["./INCAR.state_1", "./INCAR.state_2"]`。每个模板分别生成一套基准和磁扰动任务，输出中增加 `000-incar`、`001-incar` 层；即使列表只有一个文件，也会有 `000-incar`。使用单个字符串时没有这一层。

### 五种磁矩扰动模式

```json
{
  "pert_spin": [
    {"Rotation": {"angle": [45, 90], "axis": [0, 0, 1]}},
    {"Canting": {"angle": [30, 60], "seed": 12345}},
    {"Rota_Cant": {"R_angle": 30, "axis": [0, 1, 0], "C_angle": 45, "seed": 12345}},
    {"Random": {"num": 3, "seed": 12345}},
    {"Scale": {"pert": 0.10, "pert_step": 0.05}}
  ]
}
```

- `Rotation`：将每个磁矩绕给定 `axis` 按右手定则旋转 `angle` 度，保持模长。
- `Canting`：每个非零磁矩与原方向成指定 `angle` 度，方位角分别随机抽取，保持模长。
- `Rota_Cant`：固定先执行 Rotation，再执行 Canting。
- `Random`：为每个非零磁矩随机取方向并保持模长；`num` 指生成的构型数。
- `Scale`：只改变模长，按非零相对增量 `δ` 计算 `m'=(1+δ)m`；增量从 `-pert` 到 `+pert`，步长为 `pert_step`，不保留零增量。

每个 `pert_spin` 项都**从原始模板磁矩独立出发**，不会把上一项的结果作为下一项输入。输出任务数是各项变体数之和，再加一个 `origin/000000`；同一项内部的列表参数才做笛卡尔组合。`Canting`、`Rota_Cant`、`Random` 可用 `seed` 复现随机序列。旧参数 `Rcut` 和 `direction` 不再使用。

## 4. 编写 `machine.json`

沿用 `init_bulk` 的嵌套 `fp` 配置；若磁性计算要使用不同 VASP 命令，可增加同结构的 `spin` 项。以下 `remote_root`、分区、CPU 数和命令都需要按自己的集群修改：

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

未配置 `spin` 时 Stage 4 复用 `fp`。旧式扁平 `fp_*` 配置也兼容。Stage 5 的 `convert-data` 可写成示例中的单项列表，也可直接写成对象；`command` 必须为不含注释的单行命令，最后一条简单命令必须是 `nequip-data`（前面可以写 `source ... &&`，后面不要再接管道或其他命令）。程序自动在远端创建 `out/`，并在命令未指定时追加 `-p data -o out/data.extxyz`；若用户显式给出这两个参数，其值必须与这里一致。示例分区、资源和路径均需按集群修改。将 `vasp.slurm` 放入 `user_forward_files` 仅表示把文件带入任务目录；只有把命令写成 `sbatch vasp.slurm` 才会实际提交该脚本。普通 `sbatch` 入队后立即返回，可能导致 DPDispatcher 过早回传，推荐由 DPDispatcher 管理 Slurm 作业并在其中直接执行 `srun vasp_std` 或 `srun vasp_ncl`。

## 5. 运行与分阶段检查

首次使用建议逐阶段执行。每次将 `param_spin.json` 中的 `stages` 设置为对应值；Stage 4 的 `spin_action` 也按下表设置：

| 阶段 | `stages` 与命令 | 完成后检查 |
| --- | --- | --- |
| 1 结构扰动 | `[1]`；`dpgen spin_init param_spin.json` | `00.scale_pert/scale-*/000000/POSCAR` 等文件 |
| 2 AIMD | `[2]`；`dpgen spin_init param_spin.json machine.json` | 每个 `01.md` task 的 `OUTCAR`、`XDATCAR` 已回传 |
| 3 快照 | `[3]`；`dpgen spin_init param_spin.json` | `02.disp/.../00/POSCAR` 等快照 |
| 4 建任务 | `[4]` 且 `spin_action="make"`；不传 MACHINE | `03.spin/tasks.json`、`origin/000000/INCAR` 等 |
| 4 运行任务 | `[4]` 且 `spin_action="run"`；传 MACHINE | `03.spin` 每个 task 的 `OUTCAR`、`OSZICAR` |
| 5 导出数据 | `[5]`；`dpgen spin_init param_spin.json machine.json` | `04.data/selection.json`、各 scale 的 `out/data.extxyz`、raw/npy |

Stage 2 不传 MACHINE 时只建 AIMD 目录，不提交计算。Stage 4 的 `make_run` 在传 MACHINE 时建目录后提交，不传时只建目录。若机器、VASP 与 `convert-data` 已配置好，也可使用 `"stages": [1, 2, 3, 4, 5]` 和 `spin_action="make_run"` 一次执行。Stage 3 会分别检查 MD 是否完整结束、XDATCAR 是否能解析。Stage 5 不运行 VASP，但需要 MACHINE 中的 `convert-data` 配置来提交 `nequip-data`；它要求 `03.spin/tasks.json` 中所有任务已正常结束。已有的阶段输出目录不会被静默覆盖；重跑已有 Stage 4 任务请使用 `spin_action="run"`。

## 6. 输出与是否成功的判断

```text
run_spin/
├── param.json
├── 00.scale_pert/scale-1.000/000000/POSCAR
├── 01.md/scale-1.000/000000/{POSCAR,INCAR,POTCAR,OUTCAR,XDATCAR}
├── 02.disp/scale-1.000/000000/00/POSCAR
├── 03.spin/
│   ├── POTCAR
│   ├── tasks.json
│   └── scale-1.000/000000/00/origin/000000/{POSCAR,POTCAR,INCAR,OUTCAR,OSZICAR}
├── 03.spin/scale-1.000/data/{OUTCAR-1,OSZICAR-1,...}  # 相对符号链接
└── 04.data/
    ├── selection.json
    └── scale-1.000/
        ├── data -> 03.spin/scale-1.000/data
        └── out/{data.extxyz,*.raw,set/*.npy}
```

`01.md` 和 `03.spin` task 的 POSCAR/POTCAR 是真实相对 symlink，不是普通副本。`03.spin` 的磁扰动任务与 `origin/` 同级分组，例如 `000-canting/C1/`；多个 `spin_incar` 时，在 `00/` 与各组之间增加 `000-incar/` 等层。目录中的 `OUTCAR`、`XDATCAR`、`OSZICAR` 只有计算完成且回传后才存在。

新建 Stage 4 使用“编号在前”的目录名。`spin_action="run"` 仍可读取旧清单中的 `incar-000` 和 `Rotation-000` 等目录，但不会自动重命名或搬移已有任务。

Stage 5 不再进行 RMSE 筛选。它确认 `tasks.json` 中所有磁性任务正常完成后，将全部 OUTCAR/OSZICAR 按 scale 连续编号并送入转换；`selection.json` 保留来源 task、scale 和编号，`rejected` 固定为空。`convert-data` 生成 `out/data.extxyz`，随后 `out2npy` 一步生成并保留 raw 和 `set/*.npy`（`energy.npy` 为一维）。

详细编号、单/多 INCAR 目录差异、文件来源及 Stage 5 文件树见 [输出目录结构详解](SPIN_INIT_OUTPUT_STRUCTURE.md)。

## 7. 常见问题

- `dpgen spin_init` 不被识别：确认在当前环境安装的是本仓库，并运行 `python -m pip show dpgen` 检查安装位置。
- 导入时提示 `numpy.dtype size changed`：当前 NumPy 与 h5py 二进制版本不兼容；在同一个环境中安装相容版本。
- Bohrium 上传时报 `name 'oss2' is not defined`：安装 `dpdispatcher[bohrium]`。
- Stage 3 报缺少 XDATCAR：先确认 VASP 在计算节点生成了 XDATCAR，以及 DPDispatcher 是否把它回传到对应 `01.md` task。
- Stage 4 报磁矩数量不符：检查 POSCAR 原子数，以及 `MAGMOM` 和 `M_CONSTR` 的 `3 × 原子数` 个分量。
- Linux 之外出现 `WinError 1314`：这是 Windows 符号链接权限限制，不能通过把 POSCAR/POTCAR 改为 copy 来规避设计要求。

更多参数和计算检查细节见 [详细输入参考](doc/init/spin-init-usage.md)。
