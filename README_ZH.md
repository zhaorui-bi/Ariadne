<p align="center">
  <img src="docs/logo.png" alt="Ariadne logo" width="520">
</p>

<h1 align="center">Ariadne</h1>

<p align="center">
  <strong>面向珊瑚 TPS 挖掘与 cembrene 类合成酶（CeeSs）优先识别的<br>萜烯合成酶发现平台</strong>
</p>

<p align="center">
  <a href="./README.md">English</a> &nbsp;·&nbsp;
  <a href="./docs/index.md">文档</a> &nbsp;·&nbsp;
  <a href="#快速开始">快速开始</a> &nbsp;·&nbsp;
  <a href="#引用">引用</a>
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/Python-%E2%89%A53.9-0b132b?style=flat-square&logo=python&logoColor=white">
  <img alt="Version" src="https://img.shields.io/badge/Version-1.1.0-0f766e?style=flat-square">
  <img alt="License" src="https://img.shields.io/badge/License-MIT-c2410c?style=flat-square">
  <img alt="Visualization" src="https://img.shields.io/badge/Stage%204-PCA%2FLDA%20visualization-2563eb?style=flat-square">
  <img alt="ESM2" src="https://img.shields.io/badge/Optional-ESM2%20CeeSs-7c3aed?style=flat-square">
  <a href="./docs/index.md"><img alt="Docs" src="https://img.shields.io/badge/Docs-MkDocs-1d4ed8?style=flat-square"></a>
</p>

---

## 摘要

萜烯合成酶（TPS）产生了自然界中结构最多样的天然产物家族之一。珊瑚基因组中存在大量候选 TPS 蛋白，但要判断哪些候选真正对应 cembrene 类二萜产物，仍然需要可解释、可审计的筛选流程。

**Ariadne** 是一个四阶段平台，用于珊瑚 TPS 发现以及 **cembrene 类合成酶（CeeSs）** 的优先识别。当前 README 描述的主流程以候选分类、PCA/LDA 可视化和可选 ESM2 打分结束，不再把最终系统发育树作为最后一个阶段。

当前仓库不再携带生成好的 `data/` 或 HMM 数据包。请通过 `--reference-dir` 传入自己的参考 FASTA 目录；可选的 ESM2 打分使用根目录下的 `TPS.xlsx`。本文档所描述的流程不再要求用户额外生成或维护 HMM 目录。

## 主要贡献

- **一个参考目录，四个分析阶段。** 用户提供的参考 FASTA 目录支撑 discovery、filtering、classification 与 visualization。
- **可解释分类。** 候选序列被放入多类群 TPS 特征空间，并通过近邻参考完成标签转移，不只是输出一个不透明分数。
- **PCA/LDA 可视化解释。** 第四阶段聚焦 `embedding.svg`、`embedding_3d_sections.svg`、`embedding.tsv`、`embedding_variance.tsv` 与候选 cluster context 表。
- **基于蛋白质语言模型的 CeeSs 打分。** 当 `TPS.xlsx` 与可选 ESM 依赖可用时，Ariadne 会用 `P(CeeSs)` 对 coral-like 候选排序。
- **可审计输出。** 每次运行都会输出 TSV 表格和 SVG 图，从原始输入追踪到最终候选优先级。

## 方法

Ariadne 是一个顺序执行的四阶段流程。第四阶段可视化来自分类阶段生成的同一特征空间，因此输出位于 `03_classification/`，而不是单独的 tree 目录。

| 阶段 | 命令 | 输入 | 关键输出 |
|---|---|---|---|
| 1 · 发现 | `ariadne run` | 转录组 / 蛋白 FASTA | `candidates.protein.faa` |
| 2 · 过滤 | `ariadne filter` 或 `ariadne run` | 候选 FASTA | `candidates.filtered.faa` |
| 3 · 分类 | `ariadne classify` 或 `ariadne run` | 过滤后 FASTA + 参考序列 | `classification.tsv`、`nearest_neighbors.tsv` |
| 4 · PCA/LDA 可视化 | `ariadne classify` 或 `ariadne run` | 分类特征空间 | `embedding.svg`、`embedding_3d_sections.svg`、`embedding_variance.tsv` |

**第一阶段 - 发现。** Ariadne 从转录组装配或预测蛋白 FASTA 出发。转录组模式使用 Pyrodigal 预测 ORF；蛋白模式直接使用输入蛋白 FASTA。

**第二阶段 - 过滤。** 候选序列按覆盖度和最短长度过滤，并在 95% identity 下折叠近重复序列。与参考序列匹配的候选会记录到 `reference_matches.tsv`，但仍保留在 `candidates.filtered.faa` 中，避免误删 known-like 的珊瑚 TPS 等位变体。

**第三阶段 - 分类。** 过滤后的候选与参考序列被表示到一个多类群 TPS 特征空间中。Ariadne 会标准化特征矩阵，通过 *k* 近邻参考投票生成预测标签，并输出主分类表。

**第四阶段 - PCA/LDA 可视化。** 同一特征空间会被投影用于人工检查。Ariadne 优先使用有监督 LDA；当参考标签条件不足时回退到 PCA。主要结果包括 2D embedding、三视角 3D embedding、方差表与候选 cluster context。

## 安装

**推荐** - Python 3.11 + 仓库自带 Conda 环境：

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

**最小安装** - 仅核心依赖：

```bash
pip install -e .
```

**启用 ESM2 CeeSs 打分：**

```bash
pip install -e '.[esm]'   # 额外安装 torch、transformers、tqdm
```

核心依赖：`numpy >= 1.24`、`pyhmmer >= 0.12.0`、`pyrodigal >= 3.7.0`、`scikit-learn >= 1.4`、`openpyxl >= 3.1`。

## 快速开始

### 从蛋白 FASTA 运行 README 主流程

将预测蛋白 FASTA 放入例如 `my_proteins/` 的目录，然后运行：

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --skip-phylogeny \
  --output-dir results/
```

这里显式加入 `--skip-phylogeny`，让 `ariadne run` 与当前 README 的主流程保持一致：第四阶段是 PCA/LDA 可视化，而不是建树。

输出结构：

```text
results/
├── 01_discovery/          # discovery 命中与各样本蛋白 FASTA
├── 02_filtering/          # 过滤后 FASTA、filter_report.tsv、dedupe_clusters.tsv
├── 03_classification/     # 分类、近邻证据、PCA/LDA SVG 图
│                          #   （使用 TPS.xlsx 时附带 ceess_* 输出）
└── pipeline_summary.tsv
```

### 转录组模式

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --skip-phylogeny \
  --output-dir results_from_transcriptomes/
```

### 单独运行分类与可视化

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir results/03_classification/
```

### 从 Python / Notebook 调用

流程函数已在包顶层重新导出，可直接从脚本或 notebook 驱动；裸 `import ariadne` 保持轻量，重型可选依赖会在首次使用时惰性导入。

```python
import ariadne as ad

filtering = ad.filter_candidates(
    "candidate_tps.faa",
    "results/02_filtering",
    reference_dir="reference_fastas/",
)
classification = ad.classify_candidates(
    filtering["filtered_fasta"],
    reference_dir="reference_fastas/",
    output_dir="results/03_classification",
    ceess_xlsx="TPS.xlsx",
)
```

Notebook 示例见 [`tutorial/tutorial.ipynb`](tutorial/tutorial.ipynb)，命令行 smoke test 见 [`tutorial/run_example.sh`](tutorial/run_example.sh)。

## CeeSs 打分（可选）

当 `TPS.xlsx` 存在且已安装 ESM 依赖时，`ariadne run` 与 `ariadne classify` 会在特征空间分类后追加一层 ESM2 打分：

1. 从 `TPS.xlsx` 读取标注序列（`Name` / `Protein` / `Type` / `Species`）。
2. 对训练序列与 coral-like 候选计算冻结、均值池化的 ESM2 embedding。
3. 在训练 embedding 上训练轻量分类头。
4. 对每个候选打分，`P(CeeSs)` 为所有 CeeSs 正类标签概率之和。
5. 高于 `--ceess-threshold`（默认 `0.9`）的候选写入 `ceess_candidates.{tsv,fasta}`。

| `--ceess-classifier` | 分类头 |
|---|---|
| `mlp`（默认） | Torch MLP，cross-entropy + AdamW，类别加权 |
| `logreg` | Scikit-learn logistic regression + 标准化 |
| `contrastive` | Barlow Twins projection network + MLP head |

关键文件：`ceess_predictions.tsv`、`ceess_candidates.{tsv,fasta}`、`ceess_embedding.svg`、`ceess_model_metrics.tsv`。完整输出说明见 [docs/esm-type.md](./docs/esm-type.md)。

## 命令行

| 命令 | 用途 |
|---|---|
| `ariadne run` | 运行到分类与 PCA/LDA 可视化的端到端流程 |
| `ariadne filter` | 第二阶段：质量过滤与去冗余 |
| `ariadne classify` | 第三、四阶段：特征空间分类与可视化 |
| `ariadne prepare-references` | 从源数据准备干净的参考 FASTA |

> 全局参数（`--verbose`、`--log-file`）必须放在子命令之前：`ariadne --verbose run ...`

主流程参数表见 **[CLI Reference](./docs/cli-reference.md)**。

## 仓库结构

```text
Ariadne/
├── ariadne/           # 核心包
│   ├── search.py      # 第一阶段：候选发现
│   ├── filter.py      # 第二阶段：覆盖度、长度与近重复过滤
│   ├── embed.py       # 第三、四阶段：分类与 PCA/LDA 可视化
│   ├── model.py       # 可选 ESM2 CeeSs 打分
│   ├── data.py        # 参考数据管理
│   ├── utils.py       # 日志、FASTA/TSV I/O、序列工具
│   └── cli.py         # 命令行接口
├── TPS.xlsx           # 可选 CeeSs 训练标注表
├── docs/              # 文档站点
├── docs/logo.png      # README 与文档 logo
├── tutorial/          # 可运行脚本与 notebook 教程
├── environment.yml
└── pyproject.toml
```

## 文档

完整文档站点位于 [`docs/`](./docs/index.md)：

- [Getting Started](./docs/getting-started.md) - 安装与首次运行
- [Method](./docs/method.md) - 四阶段设计
- [Advanced Usage](./docs/advanced-usage.md) - 阈值调参、ESM2 分类头与可复现实验配置
- [Tutorials](./docs/tutorials.md) - 实用分析路径
- [CLI Reference](./docs/cli-reference.md) - 全部命令与参数
- [Outputs](./docs/outputs.md) - 流程产出的全部文件
- [CeeSs Classifier](./docs/esm-type.md) - ESM2 打分层
- [Citation](./docs/citation.md)

## 示例

可以用下面的命令做本地 smoke test：

```bash
bash tutorial/run_example.sh results/ my_proteins/ reference_fastas/
RUN_CEESS=1 bash tutorial/run_example.sh results/ my_proteins/ reference_fastas/
```

脚本需要传入蛋白输入目录和参考 FASTA 目录，并将流程限制在分类与 PCA/LDA 可视化。

## 开发

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
ruff check ariadne
```

## 引用

```bibtex
@software{jiang2026ariadne,
  author  = {Jiang, Zhaorui},
  title   = {Ariadne: A Coral-Centered Terpene Synthase Discovery and CeeSs Prioritization Platform},
  year    = {2026},
  url      = {https://github.com/zhaorui-bi/Ariadne},
  version = {1.1.0}
}
```

## 许可证

基于 [MIT License](./LICENSE) 发布。
