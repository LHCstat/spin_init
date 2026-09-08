# `spin_init` 使用说明

## 作用

`spin_init` 用于从一个已有 POSCAR 生成一组结构扰动后的 VASP AIMD 任务，并将每条
AIMD 轨迹的 XDATCAR 拆分成独立 POSCAR。生成的结构快照可作为后续自旋初始化或其他
构型采样步骤的输入。

```text
POSCAR → supercell/scale/perturb → AIMD → XDATCAR → POSCAR snapshots
```

当前功能不进行自旋倾斜、旋转或磁矩扰动，也不创建 `03.spin`。

## 命令

```bash
dpgen spin_init PARAM [MACHINE]
```

- `PARAM`：工作流参数文件，JSON 或 YAML。
- `MACHINE`：可选的 dpdispatcher machine 文件。stage 2 提供该文件时会自动提交
  VASP；不提供时只建立计算目录。

## 输入文件

| 文件 | 是否必需 | 作用 |
| --- | --- | --- |
| POSCAR | 是 | 初始结构；程序从这里扩胞、缩放和扰动 |
| INCAR.md | 是 | VASP AIMD 参数；`NSW` 决定预期轨迹步数 |
| POTCAR/赝势片段 | 是 | 按 POSCAR 元素顺序组成公共 POTCAR |
| spin-init.json | 是 | stages、扩胞、缩放、扰动和输入路径 |
| machine.json | 自动提交时必需 | dpdispatcher 机器、资源和 VASP 命令 |
| KPOINTS | 通常需要 | 由 machine.json 的 user forward files 传入各 task |

### POSCAR 的编写方式

使用标准 VASP POSCAR，至少应正确包含晶格、元素、各元素原子数、坐标类型和全部原子
坐标：

```text
Fe O
1.0
3.0 0.0 0.0
0.0 3.0 0.0
0.0 0.0 3.0
Fe O
1 1
Direct
0.0 0.0 0.0
0.5 0.5 0.5
```

输入结构默认就是 AIMD 前的基态/初始结构；`spin_init` 不执行预弛豫。

### INCAR.md 的编写方式

INCAR 必须配置成分子动力学而非结构优化。最短测试可以使用较小的 `NSW`：

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

材料相关参数必须由使用者确认。若 `NSW` 与 PARAM 中的 `md_nstep` 不同，程序使用
`NSW`。

### POTCAR 的准备方式

POTCAR 的元素顺序必须与 POSCAR 一致。可以传入一个已拼接文件：

```json
"potcars": ["./POTCAR"]
```

也可以按元素顺序列出多个片段：

```json
"potcars": ["./POTCAR_Fe", "./POTCAR_O"]
```

程序将这些片段顺序拼接为 `01.md/POTCAR`。赝势需由用户从有权使用的 VASP 赝势库
准备。

### spin-init.json 的编写方式

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

| 字段 | 编写要求 |
| --- | --- |
| `stages` | 使用 `1`、`2`、`3`；分别表示生成结构、建立/运行 AIMD、收集快照 |
| `from_poscar_path` | 初始 POSCAR 路径 |
| `out_dir` | 输出根目录；不会自动添加 `.spin_init` 等后缀 |
| `super_cell` | 三个正整数，例如 `[2, 2, 2]` |
| `scale` | 一个或多个正数，例如 `[0.98, 1.0, 1.02]` |
| `pert_numb` | 非负整数；总任务数为每个 scale 下 `pert_numb + 1` |
| `pert_box` | 非负晶胞扰动幅度 |
| `pert_atom` | 非负原子扰动幅度，单位 Å |
| `md_incar` | AIMD INCAR 路径 |
| `md_nstep` | 预期步数；INCAR 的 `NSW` 优先 |
| `potcars` | POTCAR 文件列表，顺序与 POSCAR 元素一致 |

### machine.json 的编写方式

`spin_init` 与 `init_bulk` 可以直接共用 machine.json。推荐嵌套格式：

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
      "batch_type": "Slurm",
      "number_node": 1,
      "cpu_per_node": 32,
      "group_size": 1,
      "queue_name": "partition"
    },
    "command": "srun vasp_std",
    "user_forward_files": ["/path/to/KPOINTS"]
  }
}
```

也兼容 `fp_machine`、`fp_resources`、`fp_command`、`fp_group_size`、
`fp_user_forward_files` 组成的旧式扁平格式。

KPOINTS、环境脚本或 `vasp.slurm` 通过 forward files 传输；实际执行哪个命令只由
`command`/`fp_command` 决定。OUTCAR 和 XDATCAR 已固定加入 backward files。

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
│       ├── POSCAR -> ../../../00.scale_pert/scale-1.000/000000/POSCAR
│       ├── INCAR  -> ../../INCAR
│       ├── POTCAR -> ../../POTCAR
│       ├── OUTCAR
│       └── XDATCAR
└── 02.disp/
    └── scale-1.000/000000/
        ├── 00/POSCAR
        ├── 01/POSCAR
        └── 02/POSCAR
```

- `00.scale_pert`：`000000` 为未扰动的缩放结构，后续目录为随机扰动结构。
- `01.md`：VASP AIMD 输入与回传结果。任务级 POSCAR、INCAR、POTCAR 是真实相对
  symbolic link。
- `02.disp`：最终快照，每个编号目录包含 XDATCAR 的一帧 POSCAR。

stage 3 会先独立检查 OUTCAR 是否完整，再使用 pymatgen `Xdatcar` 解析轨迹。所有任务
通过检查后才创建 `02.disp`。

完整的安装、集群配置、分阶段运行与排错步骤见仓库根目录
[SPIN_INIT_GUIDE.md](../../SPIN_INIT_GUIDE.md)。
