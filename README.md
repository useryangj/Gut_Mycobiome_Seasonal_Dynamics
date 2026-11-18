# ITS 真菌群落分析项目

本项目用于对 ITS 扩增子测序数据进行真菌群落分析，并结合季节、月份、宏观真菌组成、食物结构和气候变量开展统计分析。

> **论文公开地址：** 待补充  
> **数据公开地址：** 待补充  
> **项目状态：** 分析脚本整理中

## 项目内容

主要分析流程包括：

1. 原始双端 FASTQ 数据质控、引物去除和 DADA2 去噪。
2. ASV 序列表构建、嵌合体去除和 UNITE 分类注释。
3. 生成未标准化、稀释标准化和相对丰度特征表。
4. 在门、科、属、种水平汇总丰度。
5. 计算 α 多样性，并按季节和月份进行 Kruskal–Wallis 检验及 Dunn-BH 事后比较。
6. 基于 Bray-Curtis 和 Robust Aitchison 距离进行 PCoA、PERMANOVA 和 PERMDISP。
7. 绘制不同分类层级的群落组成、稀有化曲线、Venn 图、热图和 Sankey 图。
8. 筛选宏观真菌相关类群，并分析其在不同月份或季节中的丰度变化。
9. 使用月度食物结构、温度和湿度数据进行 RDA、方差分解、敏感性分析和核心属相关性分析。

## 文件说明

| 文件 | 作用 |
| --- | --- |
| [`ITS.R`](ITS.R) | ITS 数据主流程及下游群落统计、分类组成和可视化分析 |
| [`macro.R`](macro.R) | 宏观真菌相关属/种筛选、月度丰度汇总及可视化 |
| [`Monthly_RDA.R`](Monthly_RDA.R) | 月度真菌群落与食物结构、温度、湿度的 RDA 分析 |

## 分析流程

### 1. ITS 数据处理

`ITS.R` 使用 `dada2`、`ShortRead` 和 `Biostrings` 完成：

- 双端 FASTQ 文件读取；
- 无效碱基过滤；
- 使用 `cutadapt` 去除正向和反向引物；
- 质量过滤、错误模型学习和 DADA2 去噪；
- 双端序列合并；
- 去除嵌合体；
- 使用 UNITE 参考数据库进行分类注释；
- 导出 ASV、分类注释和不同分类水平的丰度表。

主要中间结果包括：

- `quality_control.txt`
- `feature_table.txt`
- `feature_table_relative.txt`
- `feature_table_unflatted.txt`
- `phylum_table.txt`
- `family_table.txt`
- `genus_table.txt`
- `species_table.txt`

### 2. α 多样性

脚本计算以下指标：

- Shannon 指数
- Simpson 指数
- Chao1
- ACE

分组变量包括：

- 季节：Spring、Summer、Autumn、Winter
- 月份：Jan–Dec

整体差异使用 Kruskal–Wallis 检验，组间比较使用 Dunn 检验，并采用 Benjamini–Hochberg 方法进行多重检验校正。

### 3. β 多样性

项目包含两类距离和分析设置：

- Bray-Curtis 距离：用于月度和季节 PERMANOVA、PERMDISP 及 PCoA；
- Robust Aitchison 距离：用于组成型数据视角下的群落差异分析。

相关结果通常包括：

- 整体 PERMANOVA；
- 两两 PERMANOVA；
- PERMDISP；
- PCoA 坐标和图形；
- BH 校正后的 `q` 值。

### 4. 宏观真菌分析

`macro.R` 根据 `macro_fungi.txt` 中的类群名单，从属水平和种水平丰度表中筛选对应类群，并输出：

- 宏观真菌属/种丰度表；
- 按月份和季节汇总的相对丰度；
- 柱状图；
- Sankey 图；
- 月度丰度流图。

### 5. 月度 RDA

`Monthly_RDA.R` 将 12 个月作为分析单位，将属水平真菌群落与以下环境或食物变量关联：

- Young leaves
- Mature leaves
- Fruits
- Flowers
- Invertebrate
- Seeds
- Roots
- Macro-fungi
- Barks
- Fallen leaves
- Others
- Temperature
- Humidity

脚本输出单因素 RDA、简化 RDA、边际效应、VIF、方差分解、留一月敏感性分析、循环移位敏感性分析以及核心属相关性热图。

默认输出目录为：

```text
Monthly_RDA_final1/
```

## 运行环境

建议使用 R 4.x，并安装脚本中使用的 R 包。核心依赖包括：

```r
install.packages(c(
  "vegan", "ggplot2", "ggpubr", "FSA", "cowplot",
  "ggvenn", "reshape2", "dplyr", "readxl", "tibble",
  "pheatmap", "ggrepel", "permute", "jsonlite",
  "htmlwidgets", "networkD3", "ggalluvial", "microeco"
))

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
BiocManager::install(c("dada2", "ShortRead", "Biostrings"))
```

此外，ITS 预处理流程需要：

- `cutadapt`
- UNITE ITS 参考数据库

## 运行前准备

当前脚本中仍保留了原始分析环境的绝对路径，例如：

```text
C:/Users/xuzy0/Desktop/xzy/split
D:/split
```

在其他计算机上运行前，需要根据实际目录修改：

1. 原始 FASTQ 文件目录；
2. cutadapt 可执行文件路径；
3. UNITE 参考数据库路径；
4. 特征表、宏观真菌表和环境数据路径；
5. 各分析结果的输出目录。

`Monthly_RDA.R` 还需要以下输入文件：

- `genus_table1.txt`
- `每月食物结构.xlsx`
- `每月温度湿度.xlsx`

## 推荐运行顺序

```text
1. ITS.R
2. macro.R
3. Monthly_RDA.R
```

建议每一步运行后检查质量控制表、输入样本数量、分组数量和输出文件，再继续下一步分析。

## 统计分析注意事项

- 当前月度 RDA 使用 12 个月作为观测单位，统计功效有限。
- 月份具有时间顺序，普通自由置换不能完全解决时间自相关问题。
- RDA、相关性和方差分解结果应视为探索性结果，不能单独解释为因果关系。
- PERMANOVA 显著时应结合 PERMDISP 结果判断组内离散度差异的影响。
- 相对丰度、稀释数据和原始计数数据不可混用；运行分析前应确认输入表类型。
- 若样本来自重复个体，应使用包含个体信息的设计约束置换，而不是完全自由置换。

## 论文信息

论文公开地址将在后续补充。建议补充以下信息：

- 论文标题；
- 作者；
- 期刊或预印本平台；
- DOI 或公开链接；
- 数据库登录号；
- 与本项目对应的图表或分析章节。

## 引用

使用本项目时，请同时引用项目所依赖的软件和数据库，包括 DADA2、UNITE、vegan、cutadapt 及其他实际使用的 R 包。
