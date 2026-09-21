# `out2npy` 简明用法

`out2npy` 一步将 `convert-data` 生成的多帧 `data.extxyz` 转成磁性数据集：

```text
data.extxyz → out2npy → output/*.raw + output/set/*.npy
```

需要 Python 3.9 或更新版本、NumPy，并在同一 Python 环境中安装本项目。
在项目根目录可执行 `python -m pip install .`。随后运行：

```bash
python -m dpgen.data.out2npy /path/to/data.extxyz /path/to/output
```

第二个参数是输出目录。目录不存在时会创建；如果目录已存在（例如 Stage 5
的 `out/` 中已有 `data.extxyz`），只要没有同名 raw 文件和 `set/`，也可以使用。
已有同名结果不会被覆盖。无需再运行第二个转换脚本。

输入可包含多帧，但原子数、元素及顺序必须一致。每帧的注释行需要
`Lattice`（9 个数）、`stress`（9 个数）、`energy`；`Properties` 需定义
`species`、`pos`（3 列）、`spin_length`、`initial_magmoms`（3 列）、
`spin_forces_vert`（3 列）和 `forces`（3 列）。缺失或不完整的帧会报错，
不会发布部分输出。

```text
output/
├── data.extxyz              # 输入位于此处时会保留
├── type_map.raw
├── type.raw
├── box.raw
├── coord.raw
├── energy.raw
├── force.raw
├── force_mag.raw
├── spin.raw
├── virial.raw
└── set/
    ├── box.npy
    ├── coord.npy
    ├── energy.npy
    ├── force.npy
    ├── force_mag.npy
    ├── spin.npy
    └── virial.npy
```

每帧对应数值 raw 文件的一行。`spin = spin_length × initial_magmoms`，
`force_mag = spin_forces_vert`，`virial = -体积 × stress`。
`.npy` 数据类型为 `float64`；`energy.npy` 的形状是 `(帧数,)`，
其他 `.npy` 的形状是 `(帧数, 对应列数)`。Stage 5 在 `convert-data`
回传 `data.extxyz` 后自动调用同一转换函数。
