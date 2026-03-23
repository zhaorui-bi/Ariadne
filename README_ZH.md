<p align="center">
  <img src="fig/ariadne_cvpr_logo.svg" alt="Ariadne logo" width="980">
</p>

<p align="center">
  <a href="./README.md"><strong>English</strong></a>
</p>

<p align="center">
  🧬 <strong>Ariadne</strong><br>
  一个面向珊瑚 TPS 挖掘与 CeeSs 优先识别的平台
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/Python-3.11%2B-0F172A?style=flat-square&logo=python&logoColor=white">
  <img alt="Workflow" src="https://img.shields.io/badge/Workflow-4%20Stages-0F766E?style=flat-square">
  <img alt="References" src="https://img.shields.io/badge/Default%20Reference-tree%2F-C2410C?style=flat-square">
  <img alt="ESM2" src="https://img.shields.io/badge/Optional-ESM2%20Type%20Model-7C3AED?style=flat-square">
  <img alt="MAFFT" src="https://img.shields.io/badge/MAFFT-required-14532D?style=flat-square">
  <img alt="IQ-TREE" src="https://img.shields.io/badge/IQ--TREE-required-1D4ED8?style=flat-square">
</p>

> Ariadne 是一个面向研究场景的平台，用于珊瑚 TPS 基因组挖掘、CeeSs 优先识别以及下游系统发育分析。  
> 当前版本采用更清晰的 `discovery -> filtering -> classification -> phylogeny` 四阶段流程，并默认以 `tree/` 作为参考来源。

## 🪸 Introduction

通过基因组挖掘来理性发现 terpene synthase（TPS），是获取新型天然产物骨架的一条重要路线。然而，发现一个新的 TPS 基因，远比准确识别“负责某个特定终产物”的 TPS 更容易，因此实验生化验证仍然是不可缺少的环节。

刺胞动物是海洋中能够产生多样萜类防御代谢物的重要类群，这意味着珊瑚基因组也可以成为 TPS 挖掘的重要资源。Ariadne 正是在这一背景下构建的：它不仅支持对珊瑚基因组进行全局 TPS 挖掘，也支持对产物特异性 synthase 的定向优先识别，这里统一称为 <strong>CeeSs</strong>。

这个平台的目标，是把计算筛选和实验验证真正接起来。在对应研究语境中，Ariadne 被用于发现 CeeSs 候选，并支撑后续异源表达验证，预测准确率达到 80%；进一步的系统发育分析也帮助解析不同 CeeSs 的演化轨迹，并为祖先酶定向进化和新型骨架的发现提供线索。

## 📢 News

- `2026-03-23`：仓库整体介绍已更新为新的 `CeeSs` 叙事，用于强调珊瑚 TPS 挖掘与 CeeSs 优先识别。
- `2026-03-21`：GitHub 默认主页切换为英文版，中文文档独立放在 [README_ZH.md](./README_ZH.md)。
- `2026-03-21`：流程收缩为更稳定的四阶段主线，`motif` 和 `benchmark` 功能已经移除。
- `2026-03-21`：`tree/` 现在同时作为 query HMM、TPS HMM library、分类背景和 phylogeny 背景的默认来源。

## ✨ 亮点

- 🌊 Tree-native：一个 `tree/` 目录驱动发现、分类和建树。
- 🧭 Feature-space 分类：候选序列会被投影到 TPS HMM 特征空间，并和参考集合做最近邻比较。
- 🤖 CeeSs-aware 分类：在判定为 `coral-like` 之后，Ariadne 可以结合 `TPS/TPS.xlsx` 和 ESM2 继续寻找 `cembrene A / cembrene B` 候选。
- 🌳 默认可建树：分类完成后，Ariadne 会直接输出 MAFFT alignment 和 IQ-TREE 系统发育树。
- 🧪 CLI 结构清晰：既支持一键 `run`，也支持分阶段的 `discover`、`filter`、`classify`、`phylogeny`。
- 📄 结果适合论文图：默认输出 `embedding.svg`、`embedding_3d_sections.svg`、局部 context tree 和最终 Newick tree。

## 🖼️ Results Preview

<p align="center">
  <img src="./docs/assets/latest_embedding.svg" alt="最新 Ariadne embedding 预览" width="100%">
</p>

<p align="center">
  <img src="./docs/assets/latest_embedding_3d_sections.svg" alt="最新 Ariadne 三视角 embedding 预览" width="100%">
</p>

<p align="center">
  <img src="./docs/assets/latest_tree.svg" alt="最新 Ariadne 系统发育树预览" width="100%">
</p>

<p align="center">
  以上结果同步自 <code>tmp_run_default_ceess/03_classification/</code> 与 <code>tmp_run_default_ceess/04_phylogeny/</code> 的最新默认流程输出，既展示了 <code>Candidate CeeSs</code> 和 <code>Candidate non-CeeSs</code> 在分类 embedding 中的区分，也展示了真实的下游 IQ-TREE 系统发育树结果。
</p>

## 🧠 整体流程

当前 Ariadne 的设计尽量保持紧凑：

1. `discovery`
   优先使用 `ariadne/hmm/query.hmm` 进行 discovery；如果内置 HMM 不可用，再回退到 `tree/` 中的参考 FASTA 自动构建 query HMM，并在蛋白 FASTA 或转录组 ORF 中搜索候选 TPS。
2. `filtering`
   去掉低 coverage、过短和近重复的候选序列。
3. `classification`
   将候选与参考一起打成 TPS HMM score matrix，做 embedding，并为每个候选给出最近参考来源预测。
4. `phylogeny`
   将过滤后的候选与参考序列合并，运行 MAFFT，并用 IQ-TREE 构建系统发育树。

## 🗂️ 仓库结构

```text
Ariadne/
├── ariadne/                # 核心包，内含 ariadne/hmm/ 默认 HMM
├── input/                  # 示例输入
├── TPS/                    # 用于 ESM type 分析的珊瑚 TPS 标注表
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
- `TPS/`
  可选的监督分析目录。`TPS/TPS.xlsx` 里保存了整理好的珊瑚 TPS 蛋白表，第二列是蛋白序列，第三列是产物 `Type` 标签。当前 `classification` 阶段会在需要时用它训练一个 ESM-based CeeSs 小模型。
- `tree/`
  当前最核心的参考目录。Ariadne 会从这里读取多物种 TPS FASTA，并在需要回退构建时生成 query HMM、TPS HMM library、分类背景和 phylogeny 背景。
- `ariadne/hmm/`
  由当前 `tree/` 数据集预先生成并内置在软件包中的默认 HMM 目录。Ariadne 现在会优先使用 `ariadne/hmm/query.hmm` 作为 discovery HMM，并使用 `ariadne/hmm/*.hmm` 作为默认 TPS HMM library。
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

### 可选安装：ESM 依赖

如果你希望 Ariadne 在 `classification` 阶段进一步识别 `cembrene A / cembrene B` 候选，请额外安装：

```bash
pip install -e '.[esm]'
```

这一步会补上 ESM2 所需的 `torch` 和 `transformers`。

### ESM model presets

当前 Ariadne 默认使用更大的 ESM2 checkpoint：

```text
facebook/esm2_t33_650M_UR50D
```

为了方便切换，CLI 也支持简写 preset：

| Preset | 实际解析到的模型 |
| --- | --- |
| `150M` | `facebook/esm2_t30_150M_UR50D` |
| `650M` | `facebook/esm2_t33_650M_UR50D` |
| `3B` | `facebook/esm2_t36_3B_UR50D` |
| `15B` | `facebook/esm2_t48_15B_UR50D` |

如果你习惯写 `504M`，Ariadne 也会接受，并在内部映射到当前默认的大模型 checkpoint。

示例：

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_classification/ \
  --ceess-model-name 150M
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

- 优先使用 `ariadne/hmm/query.hmm` 作为 discovery HMM
- 优先使用 `ariadne/hmm/*.hmm` 作为默认 TPS HMM library
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

### 在 classification 里继续寻找 CeeSs 候选

```bash
ariadne classify \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_classification/
```

如果仓库里存在 `TPS/TPS.xlsx`，并且已经安装了可选 ESM 依赖，那么 `classification` 会自动：

- 先圈出 `coral-like` 的候选
- 用 `TPS/TPS.xlsx` 训练一个小型 ESM type classifier
- 对这些 coral-like candidate 继续打分，判断它们是否更像 `cembrene A / cembrene B`
- 导出最终的 CeeSs 候选表和 fasta

重点结果文件包括：

- `classification.tsv`
- `ceess_predictions.tsv`
- `ceess_candidates.tsv`
- `ceess_candidates.fasta`
- `ceess_embedding.svg`
- `ceess_model_metrics.tsv`

如果你只想单独分析 `TPS/TPS.xlsx` 这张已标注表，也可以继续用：

```bash
ariadne esm-type \
  --xlsx TPS/TPS.xlsx \
  --output-dir esm_results/
```

这个独立命令会输出：

- `esm_embedding.svg`
- `esm_projection.tsv`
- `esm_predictions.tsv`
- `esm_confusion_matrix.tsv`
- `esm_metrics.tsv`

基于当前仓库中的 `TPS/TPS.xlsx`，默认配置下共有 `50` 条蛋白、`5` 个类型，当前 ESM baseline 的交叉验证准确率为 `0.74`。

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
- `04_phylogeny/phylogeny_preview.svg`
- `04_phylogeny/iqtree.treefile`
- `04_phylogeny/iqtree.iqtree`

这些文件对应下游手工检查、作图和生物学解释时最核心的 alignment 与 tree。

### Tutorial 4. 在 classification 之后查看 CeeSs 结果

完成 `ariadne classify` 之后，建议优先查看：

- `classification.tsv`
- `ceess_predictions.tsv`
- `ceess_candidates.tsv`
- `ceess_candidates.fasta`
- `ceess_embedding.svg`

在当前仓库自带的一个真实 candidate 集合上，Ariadne 先识别出了 `36` 条 `coral-like` 序列，再在 `--ceess-threshold 0.5` 下保留了 `12` 条高置信度 CeeSs 候选。

### Tutorial 5. 对整理好的 coral TPS 表单独跑 ESM type 模型

```bash
.venv/bin/python -m ariadne esm-type \
  --xlsx TPS/TPS.xlsx \
  --output-dir tmp_esm_results
```

建议优先查看：

- `tmp_esm_results/esm_embedding.svg`
- `tmp_esm_results/esm_metrics.tsv`
- `tmp_esm_results/esm_confusion_matrix.tsv`
- `tmp_esm_results/esm_predictions.tsv`

这条路线面向“已知 coral TPS 序列的监督式类型区分”，和主流程中的 de novo candidate discovery 是互补关系，不是替代关系。

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
| `--query-hmm` | `ariadne/hmm/query.hmm` | 显式指定 discovery HMM |
| `--tps-hmm-dir` | `ariadne/hmm/` | 显式指定 TPS HMM library |
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
| `--tps-hmm-dir` | `ariadne/hmm/` | TPS HMM library 目录 |
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
│   ├── phylogeny_preview.svg
│   ├── phylogeny_sequence_map.tsv
│   ├── iqtree.treefile
│   └── iqtree.iqtree
└── pipeline_summary.tsv
```

## 🔬 建议的结果阅读顺序

1. 先看 `pipeline_summary.tsv`，确认流程是否完整结束。
2. 再看 `03_classification/classification.tsv`，理解每条 candidate 的预测来源。
3. 打开 `03_classification/embedding.svg`，检查全局分布。
4. 建议先看 `04_phylogeny/phylogeny_preview.svg`，再结合 `04_phylogeny/iqtree.treefile` 和 `04_phylogeny/iqtree.iqtree` 做完整系统发育解释。

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
