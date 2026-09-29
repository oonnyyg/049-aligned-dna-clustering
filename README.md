# Aligned DNA Clustering

本地相似性分析工具：读取一批**已经对齐**的等长 DNA 序列（FASTA），
产出可复核的逐对比较表、对称距离矩阵和 UPGMA 层次聚类树。

- 仅使用 Python 3.12 标准库，不调用任何在线服务或第三方计算库。
- 不重新比对、不推断真实进化关系；`N` 与间隙（`-`）表示缺失/不可靠位点，
  其影响对每一对序列单独计算。

## 用法

```powershell
python -m alignment_cluster examples\demo.fasta -o examples\demo_output
```

参数：

- `fasta`：输入的已对齐 FASTA 文件。
- `-o/--output-dir`：报告输出目录（不存在会创建）。

成功时在输出目录生成三个文件；输入非法时进程以非零状态退出，
**不会生成或覆盖任何报告文件**（所有校验在写盘之前完成）。

## 输入规则

- 至少两条序列，全部非空且长度一致。
- 标题行以 `>` 开始；其后第一个非空白词为标识，必须非空且全文件唯一，
  其余文字作为描述（可为空）。
- 序列允许跨行，文件中允许空白行；序列行内空白被忽略，字母统一转大写。
- 只接受 `A`、`C`、`G`、`T`、`N`、`-`。
- 以下情况一律拒绝：标题之前出现序列、重复标识、空标识、
  非法字符、序列长度不一致、序列为空、序列少于两条。

## 计算口径

### 逐对距离（缺失位点的影响）

对每一个无序序列对，**只统计双方都为 A/C/G/T 的列**：

- 某列任一方为 `N` 或 `-`，该列仅从**当前这一对**中排除，
  不会影响其他序列对，也不会从任何序列中删除位点。
- `comparable_sites` = 双方均为 A/C/G/T 的列数。
- `differences` = 其中碱基不同的列数。
- `distance = differences / comparable_sites`。

若任一对的 `comparable_sites` 为 0，整次分析失败，
  错误信息列出对应的两个标识（可能同时列出多对）。

### UPGMA 聚类

- 初始每个样本为一个单序列簇，保留全部叶节点（零距离样本不会被合并删除）。
- 每轮选择簇间距离最小的两簇合并；簇间距离是跨两簇
  **全部原始叶节点对距离的算术平均**（等价于按簇大小加权，
  大小不同的簇不会被等权处理）。
- 距离相等时，按两簇各自最小标识组成的有序对取字典序；
  子树也按各自最小标识从左到右排列。
- 叶节点高度为 0；合并节点高度 = 本次簇间距离 / 2；
  分支长度 = 父高度 − 子高度（非负，根的分支长度为 0）。
- 全部比较与决策使用精确分数（`fractions.Fraction`），
  **不会用提前舍入的小数决定合并顺序**；`value` 字段的浮点值仅为便于阅读。
- 输出不依赖输入记录顺序：重排 FASTA 得到相同的矩阵与树。

## 输出格式

### `distance_matrix.csv`

- 行列均按标识字典序排列；矩阵对称、对角线为 0。

### `pairwise_comparisons.csv`

- 每个无序序列对仅一行，按 `(id_a, id_b)` 字典序排列。
- 列：`id_a, id_b, comparable_sites, differences, distance`。
- `distance` 为浮点展示，精确值见 `cluster_tree.json` 中的
  `numerator/denominator`。

### `cluster_tree.json`

- `identifiers`：字典序标识列表。
- `pairwise_distances`：每个 `"id_a|id_b"` 的精确分数与浮点值。
- `merges`：按发生顺序记录每次合并，含左右两簇的全部成员
  （`members`、`min_identifier`）、精确簇间距离 `cluster_distance`
  和节点高度 `node_height`，可独立还原每次合并。
- `tree`：最终嵌套树，叶节点 `type: "leaf"`，合并节点
  `type: "internal"`，每个节点带精确 `height` 与 `branch_length`。

## 测试

```powershell
python -m unittest discover -s tests
```

覆盖：缺失位点逐对掩码、无可比较位点失败、零距离样本保留、
并列合并裁决、非等大小簇的平均距离、输入换序不变性、
分支长度与高度关系，以及非法输入非零退出且不覆盖已有报告。
