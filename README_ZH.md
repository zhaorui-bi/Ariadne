<p align="center">
  <img src="fig/ariadne_cvpr_logo.svg" alt="Ariadne logo" width="980">
</p>

<p align="center">
  <a href="./README.md"><strong>English</strong></a>
</p>

<p align="center">
  🧬 <strong>Ariadne</strong><br>
  一个面向珊瑚及跨物种 TPS 参考库的 tree-native terpene synthase 挖掘与系统发育分析平台
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/Python-3.11%2B-0F172A?style=flat-square&logo=python&logoColor=white">
  <img alt="Workflow" src="https://img.shields.io/badge/Workflow-4%20Stages-0F766E?style=flat-square">
  <img alt="References" src="https://img.shields.io/badge/Default%20Reference-tree%2F-C2410C?style=flat-square">
  <img alt="MAFFT" src="https://img.shields.io/badge/MAFFT-required-14532D?style=flat-square">
  <img alt="IQ-TREE" src="https://img.shields.io/badge/IQ--TREE-required-1D4ED8?style=flat-square">
</p>

> Ariadne 是一个面向研究场景的命令行平台，用于 terpene synthase 候选发现、TPS 特征空间分类、多序列比对与系统发育重建。  
> 当前版本采用更清晰的 `discovery -> filtering -> classification -> phylogeny` 四阶段流程，并默认以 `tree/` 作为参考来源。

## 📢 News

- `2026-03-21`：GitHub 默认主页切换为英文版，中文文档独立放在 [README_ZH.md](./README_ZH.md)。
- `2026-03-21`：流程收缩为更稳定的四阶段主线，`motif` 和 `benchmark` 功能已经移除。
- `2026-03-21`：`tree/` 现在同时作为 query HMM、TPS HMM library、分类背景和 phylogeny 背景的默认来源。

## ✨ 亮点

- 🌊 Tree-native：一个 `tree/` 目录驱动发现、分类和建树。
- 🧭 Feature-space 分类：候选序列会被投影到 TPS HMM 特征空间，并和参考集合做最近邻比较。
- 🌳 默认可建树：分类完成后，Ariadne 会直接输出 MAFFT alignment 和 IQ-TREE 系统发育树。
- 🧪 CLI 结构清晰：既支持一键 `run`，也支持分阶段的 `discover`、`filter`、`classify`、`phylogeny`。
- 📄 结果适合论文图：默认输出 `embedding.svg`、`embedding_3d_sections.svg`、局部 context tree 和最终 Newick tree。

## 🧠 整体流程

当前 Ariadne 的设计尽量保持紧凑：

1. `discovery`
   从 `tree/` 中的参考 FASTA 自动构建 query HMM，并在蛋白 FASTA 或转录组 ORF 中搜索候选 TPS。
2. `filtering`
   去掉低 coverage、过短和近重复的候选序列。
3. `classification`
   将候选与参考一起打成 TPS HMM score matrix，做 embedding，并为每个候选给出最近参考来源预测。
4. `phylogeny`
   将过滤后的候选与参考序列合并，运行 MAFFT，并用 IQ-TREE 构建系统发育树。

## 🗂️ 仓库结构

```text
Ariadne/
├── ariadne/                # 核心包
├── input/                  # 示例输入
├── tree/                   # 默认参考 FASTA 集合
├── output/                 # 历史输出示例，仅保留作参考
├── fig/                    # logo 和 figures
├── environment.yml         # conda 环境
├── install.sh
└── pyproject.toml
```

## 📁 数据约定

- `input/`
  默认示例输入目录，通常作为 `--protein-folder` 的来源。
- `tree/`
  当前最核心的参考目录。Ariadne 会从这里读取多物种 TPS FASTA，并自动构建 query HMM、TPS HMM library、分类背景和 phylogeny 背景。
- `output/`
  历史示例输出目录，只用于保留旧结果展示，不再参与当前默认流程。

## 🛠️ 安装

### Conda

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

### venv

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
bash install.sh
```

### 手动安装

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

## 📦 依赖

- Python `>= 3.9`
- `mafft`
- `iqtree` 或 `iqtree2`
- `numpy >= 1.24`
- `openpyxl >= 3.1`
- `pyhmmer >= 0.12.0`
- `pyrodigal >= 3.7.0`
- `scikit-learn >= 1.4`

推荐使用 [environment.yml](./environment.yml) 中的 Python `3.11` 环境。

## 🚀 快速开始

### 一键完整流程

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

这条命令会自动完成：

- 从 `tree/` 自动构建 discovery query HMM
- 从 `tree/` 自动构建 TPS HMM library
- 发现候选蛋白
- 过滤和去冗余
- 在 TPS 特征空间中完成分类
- 直接生成 alignment 和 phylogeny

### 转录组模式

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results_from_transcriptomes/
```

### 只做分类

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_classification/
```

### 只做比对和建树

```bash
ariadne phylogeny \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_phylogeny/
```

## 🎓 Tutorial

### Tutorial 1. 直接跑仓库自带示例

```bash
.venv/bin/python -m ariadne run \
  --protein-folder input \
  --reference-dir tree \
  --output-dir tmp_run_example
```

预期会得到：

- `tmp_run_example/01_discovery/`
- `tmp_run_example/02_filtering/`
- `tmp_run_example/03_classification/`
- `tmp_run_example/04_phylogeny/`
- `tmp_run_example/pipeline_summary.tsv`

### Tutorial 2. 先看分类结果

优先打开：

- `03_classification/classification.tsv`
- `03_classification/nearest_neighbors.tsv`
- `03_classification/embedding.svg`
- `03_classification/embedding_3d_sections.svg`

这些文件可以帮助你判断每个 candidate 最接近哪个参考来源、置信度如何，以及它在 TPS 特征空间里落在哪里。

### Tutorial 3. 再看系统发育结果

优先打开：

- `04_phylogeny/phylogeny_input.fasta`
- `04_phylogeny/phylogeny_alignment.fasta`
- `04_phylogeny/iqtree.treefile`
- `04_phylogeny/iqtree.iqtree`

这些文件对应下游手工检查、作图和生物学解释时最核心的 alignment 与 tree。

## 🧪 CLI Reference

### 主命令

| 命令 | 作用 |
| --- | --- |
| `ariadne run` | 一键完整流程 |
| `ariadne discover` | 从蛋白或转录组做 HMM 候选发现 |
| `ariadne filter` | coverage、长度和近重复过滤 |
| `ariadne classify` | TPS 特征空间分类 |
| `ariadne phylogeny` | MAFFT 比对 + IQ-TREE 建树 |
| `ariadne build-hmm` | 从一个 FASTA/MSA 源构建单个 HMM |
| `ariadne build-tps-hmm-library` | 从多个参考 FASTA/MSA 构建 TPS HMM library |

### `ariadne run`

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--protein-folder` | `None` | 蛋白 FASTA 目录 |
| `--transcriptomes` | `None` | 转录组 FASTA 输入 |
| `--reference-dir` | 必填 | tree-native 参考目录 |
| `--output-dir` | 必填 | 输出根目录 |
| `--query-hmm` | 自动构建 | 显式指定 discovery HMM |
| `--tps-hmm-dir` | 自动构建 | 显式指定 TPS HMM library |
| `--min-coverage` | `10.0` | 过滤 coverage 阈值 |
| `--min-length` | `300` | 最小蛋白长度 |
| `--identity-threshold` | `0.95` | 近重复合并阈值 |
| `--top-k` | `5` | 最近邻投票数 |
| `--tree-neighbors` | `12` | 每个局部 context tree 使用的参考邻居数 |
| `--skip-phylogeny` | `False` | 跳过 MAFFT + IQ-TREE |
| `--mafft-mode` | `--auto` | MAFFT 模式 |
| `--iqtree-model` | `LG` | IQ-TREE 模型 |
| `--iqtree-threads` | `AUTO` | IQ-TREE 线程 |
| `--iqtree-bootstrap` | `None` | 可选 bootstrap |
| `--no-iqtree-fast` | `False` | 关闭默认 `--fast` |

### `ariadne discover`

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--protein-folder` | `None` | 蛋白 FASTA 目录 |
| `--transcriptomes` | `None` | 转录组 FASTA 输入 |
| `--protein-glob` | 常见 FASTA 模式 | 在 `--protein-folder` 下递归寻找文件 |
| `--hmm` | 必填 | discovery HMM |
| `--output-dir` | 必填 | discovery 输出目录 |
| `--min-score` | `None` | 最小 HMM score |
| `--max-evalue` | `None` | 最大 E-value |

### `ariadne filter`

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--input-fasta` | 必填 | 候选蛋白 FASTA |
| `--output-dir` | 必填 | filter 输出目录 |
| `--min-coverage` | `10.0` | 去掉低 coverage 候选 |
| `--min-length` | `300` | 去掉短蛋白 |
| `--identity-threshold` | `0.95` | 近重复阈值 |

### `ariadne classify`

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--candidates` | 必填 | 过滤后的 candidate FASTA |
| `--reference-dir` | 必填 | 参考 FASTA 目录 |
| `--output-dir` | 必填 | classification 输出目录 |
| `--tps-hmm-dir` | 自动构建 | TPS HMM library 目录 |
| `--top-k` | `5` | 最近邻投票数 |
| `--tree-neighbors` | `12` | 局部 context tree 使用的邻居数 |

### `ariadne phylogeny`

| 参数 | 默认值 | 说明 |
| --- | --- | --- |
| `--candidates` | 必填 | 过滤后的 candidate FASTA |
| `--reference-dir` | 必填 | 参考 FASTA 目录 |
| `--output-dir` | 必填 | phylogeny 输出目录 |
| `--mafft-bin` | 自动识别 | 显式指定 MAFFT 可执行文件 |
| `--mafft-mode` | `--auto` | MAFFT 模式 |
| `--iqtree-bin` | 自动识别 | 显式指定 IQ-TREE 可执行文件 |
| `--iqtree-model` | `LG` | 替换模型 |
| `--iqtree-threads` | `AUTO` | 线程设置 |
| `--iqtree-bootstrap` | `None` | 可选 bootstrap |
| `--no-iqtree-fast` | `False` | 关闭默认 fast 模式 |

## 📤 输出结构

```text
results/
├── 01_discovery/
│   ├── candidates.protein.faa
│   ├── candidates.orf.fna
│   └── candidates.hits.tsv
├── 02_filtering/
│   ├── candidates.filtered.faa
│   ├── filter_report.tsv
│   ├── dedupe_clusters.tsv
│   └── manual_review.tsv
├── 03_classification/
│   ├── tps_features.tsv
│   ├── embedding.tsv
│   ├── embedding.svg
│   ├── embedding_3d_sections.svg
│   ├── classification.tsv
│   ├── nearest_neighbors.tsv
│   ├── candidate_cluster_context.tsv
│   ├── global_context_tree.nwk
│   └── trees/
├── 04_phylogeny/
│   ├── phylogeny_input.fasta
│   ├── phylogeny_alignment.fasta
│   ├── phylogeny_sequence_map.tsv
│   ├── iqtree.treefile
│   └── iqtree.iqtree
└── pipeline_summary.tsv
```

## 🔬 建议的结果阅读顺序

1. 先看 `pipeline_summary.tsv`，确认流程是否完整结束。
2. 再看 `03_classification/classification.tsv`，理解每条 candidate 的预测来源。
3. 打开 `03_classification/embedding.svg`，检查全局分布。
4. 最后结合 `04_phylogeny/iqtree.treefile` 和 `04_phylogeny/iqtree.iqtree` 做系统发育解释。

## ⚠️ 说明

- `tree/` 现在是唯一的默认参考入口。
- `Alignment.fasta` 风格的旧入口已经不属于当前主流程。
- `motif` 与 `benchmark` 功能是有意移除的，以保持当前版本聚焦和稳定。
- 仓库里仍然保留 `output/` 历史示例，但当前流程不依赖它。

## 📚 引用

如果 Ariadne 对你的工作有帮助，建议至少引用下面的软件条目；之后如果你发布了 DOI 或正式论文，也可以再替换成对应版本。

```bibtex
@software{jiang2026ariadne,
  author       = {Jiang, Zhaorui},
  title        = {Ariadne: A Tree-Native Terpene Synthase Discovery and Phylogeny Platform},
  year         = {2026},
  url          = {https://github.com/zhaoruijiang26/Ariadne},
  version      = {0.1.0}
}
```

## 🤝 致谢

Ariadne 的目标是把 TPS candidate mining 和下游的 phylogenetic interpretation 连接起来，特别适合需要跨物种参考背景的 coral-centered discovery 场景。
