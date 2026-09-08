# `spin_init` 中文上手手册

本文介绍如何安装并运行 DP-GEN 的 `spin_init` 工作流。更紧凑的功能与输入输出说明见
[spin_init 使用说明](doc/init/spin-init-usage.md)。

## 1. 功能与适用范围

`spin_init` 从一个已有 VASP POSCAR 出发，生成结构扰动后的短程 AIMD 任务，并把
XDATCAR 中的每一帧导出成独立 POSCAR：

```text
初始 POSCAR
  → 扩胞、缩放、晶胞扰动和原子扰动
  → 多个 POSCAR
  → VASP AIMD
  → OUTCAR + XDATCAR 回传
  → 逐帧 POSCAR
```

当前版本只负责结构扰动、AIMD 和轨迹快照导出。它尚不生成自旋倾斜、自旋旋转或
`03.spin` 计算任务；工作流名称表示这些快照将作为后续自旋初始化的结构输入。

`spin_init` 是独立命令，不会改变原来的 `dpgen init_bulk`。

## 2. 安装

建议在 Linux 计算集群的 Python 3.9 或更高版本环境中使用：

```bash
git clone https://github.com/LHCstat/spin_init.git
cd spin_init

python3.9 -m venv venv
source venv/bin/activate
python -m pip install --upgrade pip
pip install -e ".[test]"
```

验证命令是否已经注册：

```bash
dpgen spin_init -h
```

运行单元测试：

```bash
python -m unittest tests.data.test_spin_init -v
```

正式任务要求 Linux 支持 symbolic link。Windows 可以运行不依赖真实符号链接的纯
Python 测试，但不应把任务级 POSCAR/POTCAR 改成普通副本来规避权限问题。

## 3. 命令格式

```bash
dpgen spin_init PARAM [MACHINE]
```

- `PARAM`：必需，JSON 或 YAML 工作流参数。
- `MACHINE`：可选，DP-GEN/dpdispatcher 的机器与提交参数。

没有 `MACHINE` 时，stage 2 只建立 AIMD 目录，不提交任务。提供 `MACHINE` 且运行
stage 2 时，程序会在建好目录后立即调用 dpdispatcher 提交 VASP。

## 4. 需要准备的文件

建议在运行目录中准备：

```text
work/
├── POSCAR
├── INCAR.md
├── POTCAR
├── KPOINTS              # 通常需要；通过 machine.json 转发
├── spin-init.json
└── machine.json         # 需要自动提交 VASP 时使用
```

### 4.1 POSCAR

POSCAR 是扩胞、缩放和扰动的起点，可以是原胞、常规胞或已优化结构。元素名称、元素
顺序和原子数必须正确。例如：

```text
Fe O
1.0
3.000000 0.000000 0.000000
0.000000 3.000000 0.000000
0.000000 0.000000 3.000000
Fe O
1 1
Direct
0.000000 0.000000 0.000000
0.500000 0.500000 0.500000
```

程序不会先执行结构弛豫。请自行确认输入结构适合作为 AIMD 初始结构。

### 4.2 INCAR.md

该文件必须是 VASP AIMD 的 INCAR，而不是离子弛豫 INCAR。以下仅是短测试模板，实际
参数应根据材料、赝势和计算规范调整：

```text
SYSTEM = spin_init_test
PREC   = Normal
ENCUT  = 500
EDIFF  = 1E-5

IBRION = 0
NSW    = 3
POTIM  = 1.0
TEBEG  = 300
TEEND  = 300
SMASS  = 0
ISIF   = 2

ISMEAR = 1
SIGMA  = 0.1
LWAVE  = F
LCHARG = F
```

注意：

- `NSW` 是 AIMD 步数。
- 如果 `NSW` 与 `spin-init.json` 中的 `md_nstep` 不同，程序沿用 `init_bulk`
  的规则，以 `NSW` 为准并打印 warning。
- 极小集成测试可设置 `NSW = 2` 或 `3`。
- `IBRION`、恒温器、时间步长和电子收敛参数应符合所用 VASP 版本及体系要求。

### 4.3 POTCAR

POTCAR 必须来自用户有权使用的 VASP 赝势库，本仓库不提供赝势内容。

`potcars` 可以写一个已拼接 POTCAR，也可以按 POSCAR 的元素顺序列出多个片段：

```json
"potcars": ["./POTCAR_Fe", "./POTCAR_O"]
```

程序会按列表顺序拼接成公共 `01.md/POTCAR`。因此列表顺序必须与 POSCAR 的元素顺序
一致。

### 4.4 KPOINTS 和其他 VASP 输入

`spin_init` 不自动生成 KPOINTS。可以准备一个适合超胞 AIMD 的 KPOINTS，例如仅用于
快速测试的 Gamma 点：

```text
Automatic mesh
0
Gamma
1 1 1
0 0 0
```

通过 machine.json 的 `user_forward_files` 或 `fp_user_forward_files` 将 KPOINTS 及
其他附加文件放入每个 AIMD task。

## 5. 编写 spin-init.json

最小示例：

```json
{
  "stages": [1, 2, 3],
  "from_poscar_path": "./POSCAR",
  "out_dir": "./run_spin",
  "super_cell": [1, 1, 1],
  "scale": [1.0],
  "pert_numb": 1,
  "pert_box": 0.03,
  "pert_atom": 0.01,
  "md_incar": "./INCAR.md",
  "md_nstep": 3,
  "potcars": ["./POTCAR"]
}
```

参数说明：

| 参数 | 类型 | 说明 |
| --- | --- | --- |
| `stages` | 整数列表 | `1` 生成结构，`2` 建立/运行 AIMD，`3` 导出快照 |
| `from_poscar_path` | 字符串 | 初始 POSCAR 路径 |
| `out_dir` | 字符串 | 输出根目录；默认 `.`，不会自动添加后缀 |
| `super_cell` | 3 个正整数 | 三个晶格方向的扩胞倍数 |
| `scale` | 正数列表 | 各向同性晶格缩放值，每个值对应一个 `scale-*` 目录 |
| `pert_numb` | 非负整数 | 每个 scale 的随机扰动数，不包含 `000000` |
| `pert_box` | 非负浮点数 | 晶胞扰动幅度，与 `init_bulk` 算法一致 |
| `pert_atom` | 非负浮点数 | 原子坐标扰动上限，单位 Å |
| `md_incar` | 字符串 | AIMD INCAR 路径 |
| `md_nstep` | 非负整数 | 预期 AIMD 步数；INCAR 中的 `NSW` 优先 |
| `potcars` | 字符串列表 | POTCAR 或 POTCAR 片段，按元素顺序排列 |

每个 scale 会生成 `pert_numb + 1` 个结构：

- `000000`：只扩胞和缩放，不施加随机扰动。
- `000001` 到 `00000N`：随机晶胞/原子扰动结构。

## 6. 编写 machine.json

`spin_init` 和 `init_bulk` 可以使用同一份 machine.json。推荐使用当前嵌套格式，也兼容
旧式扁平 `fp_*` 格式。

### 6.1 推荐的嵌套格式

```json
{
  "api_version": "1.0",
  "fp": {
    "machine": {
      "batch_type": "Slurm",
      "context_type": "local",
      "local_root": "./",
      "remote_root": "/path/to/your/dpdispatcher/work"
    },
    "resources": {
      "batch_type": "Slurm",
      "number_node": 1,
      "cpu_per_node": 32,
      "group_size": 1,
      "queue_name": "your_partition",
      "module_list": ["your_vasp_module"],
      "source_list": ["/path/to/vasp_env.sh"],
      "envs": {},
      "custom_flags": ["#SBATCH -t 00:10:00"]
    },
    "command": "srun vasp_std",
    "user_forward_files": ["/path/to/KPOINTS"],
    "user_backward_files": []
  }
}
```

### 6.2 与旧式 init_bulk 相同的扁平格式

```json
{
  "api_version": "1.0",
  "fp_machine": {
    "batch_type": "Slurm",
    "context_type": "local",
    "local_root": "./",
    "remote_root": "/path/to/your/dpdispatcher/work"
  },
  "fp_resources": {
    "batch_type": "Slurm",
    "number_node": 1,
    "cpu_per_node": 32,
    "queue_name": "your_partition",
    "module_list": ["your_vasp_module"],
    "source_list": ["/path/to/vasp_env.sh"]
  },
  "fp_command": "srun vasp_std",
  "fp_group_size": 1,
  "fp_user_forward_files": ["/path/to/KPOINTS"],
  "fp_user_backward_files": []
}
```

填写时注意：

- `local_root` 使用 `"./"`。
- `remote_root` 必须是当前用户可读写的实际目录。
- 队列字段通常是 `queue_name`；具体字段以安装的 dpdispatcher 版本为准。
- `module_list` 中填写 module 名称；`source_list` 中填写需要 `source` 的脚本路径。
- `command`/`fp_command` 决定实际运行的命令。
- 把 `vasp.slurm` 放入 forward files 只会把文件带入 task，不代表自动执行它。只有
  command 明确写成 `sbatch vasp.slurm` 时才会执行该命令。
- 使用 dpdispatcher 的 Slurm backend 时，一般让 dpdispatcher 生成并提交作业，command
  填作业内执行的 VASP 命令，例如 `srun vasp_std`。

`OUTCAR` 和 `XDATCAR` 已由 `spin_init` 强制加入 backward files，无需用户重复填写。
其他需要回传的文件可放入用户 backward files。

## 7. 三阶段运行方式

### 7.1 只生成扰动结构

将参数设置为：

```json
"stages": [1]
```

运行：

```bash
dpgen spin_init spin-init.json
```

生成 `00.scale_pert`，不需要 VASP 或 machine.json。

### 7.2 只建立 AIMD 目录，不提交

在 stage 1 已完成的前提下，将参数设置为：

```json
"stages": [2]
```

不传 MACHINE：

```bash
dpgen spin_init spin-init.json
```

程序建立 `01.md` 和真实相对符号链接，但不会提交 VASP。

### 7.3 建立并自动提交 AIMD

```bash
dpgen spin_init spin-init.json machine.json
```

只要 `stages` 中包含 `2` 且提供了 MACHINE，程序就会通过 dpdispatcher 提交所有 AIMD
tasks。任务结束后，dpdispatcher 将 OUTCAR 和 XDATCAR 回传到各任务目录。

### 7.4 只收集已经完成的轨迹

把 `stages` 改为：

```json
"stages": [3]
```

运行：

```bash
dpgen spin_init spin-init.json
```

stage 3 不提交任务，只检查本地 `01.md` 中的 OUTCAR/XDATCAR 并生成 `02.disp`。即使只
运行 stage 3，也应保留可读取的 `md_incar`，以确定 `NSW`。

### 7.5 一次运行完整流程

```json
"stages": [1, 2, 3]
```

```bash
dpgen spin_init spin-init.json machine.json
```

适合 machine 配置已经验证通过的场景。首次部署建议按 stage 1、stage 2、stage 3 分开
检查。

## 8. 输出目录

示例输出：

```text
run_spin/
├── param.json
├── 00.scale_pert/
│   └── scale-1.000/
│       ├── 000000/POSCAR
│       └── 000001/POSCAR
├── 01.md/
│   ├── INCAR
│   ├── POTCAR
│   └── scale-1.000/
│       ├── 000000/
│       │   ├── POSCAR -> ../../../00.scale_pert/scale-1.000/000000/POSCAR
│       │   ├── INCAR  -> ../../INCAR
│       │   ├── POTCAR -> ../../POTCAR
│       │   ├── OUTCAR
│       │   └── XDATCAR
│       └── 000001/
└── 02.disp/
    └── scale-1.000/
        ├── 000000/
        │   ├── 00/POSCAR
        │   ├── 01/POSCAR
        │   └── 02/POSCAR
        └── 000001/
```

三个阶段的含义：

- `00.scale_pert`：扩胞、缩放和结构扰动结果。
- `01.md`：每个扰动结构对应的 VASP AIMD task。POSCAR、INCAR、POTCAR 为相对符号
  链接；OUTCAR/XDATCAR 在计算完成并回传后出现。
- `02.disp`：最终轨迹快照。每个 `NN/POSCAR` 对应 XDATCAR 中的一帧。

工作流不会静默覆盖已有阶段目录。如需重新生成某个阶段，应先备份并明确处理对应目录。

## 9. 完成性与解析检查

stage 3 分别执行两类检查：

1. **AIMD 完成性**：只检查 OUTCAR，要求存在正常结束的 `Elapse` 标志，并要求
   `TOTAL-FORCE` 块数与 `NSW`/`md_nstep` 一致。
2. **轨迹可解析性**：只检查 XDATCAR，要求文件存在、非空、可由 pymatgen `Xdatcar`
   解析，并要求各帧原子数和元素顺序一致。

所有 task 预检查通过后才创建 `02.disp`，避免部分成功造成不完整数据集。

## 10. 最小真实测试建议

首次在超算验证时建议：

```json
{
  "scale": [1.0],
  "pert_numb": 1,
  "md_nstep": 3
}
```

同时在 INCAR 中设置 `NSW = 3`。依次确认：

1. `00.scale_pert` 中有 `000000` 和 `000001`。
2. `01.md` 中 POSCAR、INCAR、POTCAR 是真实相对符号链接。
3. VASP 计算节点产生 OUTCAR 和 XDATCAR。
4. dpdispatcher 将两个文件回传到本地 task。
5. stage 3 在 `02.disp` 下为所有 XDATCAR 帧生成 POSCAR。

## 11. 常见问题

### `02.disp` 没有生成

查看错误指出的 scale、perturbation task 和文件。常见原因包括 OUTCAR 未正常结束、实际
完成步数不足、XDATCAR 未回传或轨迹文件损坏。

### task 中没有 KPOINTS

把 KPOINTS 加到 machine.json 的 user forward files。`spin_init` 不自动生成 KPOINTS。

### `md_nstep` 与实际步数不一致

检查 INCAR 的 `NSW`。程序以 `NSW` 为准。

### Windows 报 symbolic link 权限错误

正式工作流应在支持真实符号链接的 Linux 环境运行。不要把 symlink 改为 copy，也不要
用这个环境错误判断结构生成或 XDATCAR parser 失败。

### 转发了 vasp.slurm 但没有执行

forward files 只负责传输文件。实际执行内容由 machine.json 的 command 决定。

### 想重新运行某个阶段

阶段目录已存在时程序会拒绝覆盖。先备份需要保留的结果，再明确移动或删除对应的
`00.scale_pert`、`01.md` 或 `02.disp`。

## 12. 快速验收清单

- [ ] `dpgen spin_init -h` 能显示帮助。
- [ ] stage 1 生成 `pert_numb + 1` 个 POSCAR。
- [ ] stage 2 的 POSCAR、INCAR、POTCAR 是真实相对符号链接。
- [ ] machine.json 可同时用于 `spin_init` 和 `init_bulk`。
- [ ] 每个 AIMD task 回传 OUTCAR 和 XDATCAR。
- [ ] stage 3 为 XDATCAR 的每一帧生成独立 POSCAR。
- [ ] 当前流程没有生成 `03.spin`，也没有施加自旋扰动。
