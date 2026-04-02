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
  <img alt="Python" src="https://img.shields.io/badge/Python-%3E%3D3.9-0F172A?style=flat-square&logo=python&logoColor=white">
  <img alt="Workflow" src="https://img.shields.io/badge/Workflow-4%20Stages-0F766E?style=flat-square">
  <img alt="Reference Root" src="https://img.shields.io/badge/Reference-tree%2F-C2410C?style=flat-square">
  <img alt="CeeSs" src="https://img.shields.io/badge/Optional-ESM2%20CeeSs-7C3AED?style=flat-square">
  <img alt="MAFFT" src="https://img.shields.io/badge/MAFFT-required%20for%20phylogeny-14532D?style=flat-square">
  <img alt="IQ-TREE" src="https://img.shields.io/badge/IQ--TREE-required%20for%20phylogeny-1D4ED8?style=flat-square">
  <a href="./docs/index.md"><img alt="Documentation" src="https://img.shields.io/badge/Docs-Project%20Documentation-2563EB?style=flat-square&logo=readthedocs&logoColor=white"></a>
</p>

> Ariadne 是一个面向研究场景的平台，用于珊瑚 TPS 挖掘、CeeSs 优先识别以及下游系统发育分析。  
> 当前版本采用 tree-native 的四阶段主线流程：`discovery -> filtering -> classification -> phylogeny`。

## 简介

发现一个新的 terpene synthase 很重要，但更难的问题通常是识别“哪个 synthase 对应某个特定产物”。Ariadne 就是围绕这个问题构建的。

这个平台聚焦于珊瑚 TPS 发现和 CeeSs 优先识别，结合了 HMM 引导的候选筛选、TPS 特征空间分类、可选的 ESM 辅助 `cembrene A / cembrene B` 候选打分，以及直接面向 MAFFT + IQ-TREE 的系统发育分析流程。

它的核心设计原则很简单：同一个经过整理的 `tree/` 参考目录，应该同时服务于 discovery、候选解释和最终系统发育定位，让整个流程始终处在一致的生物学背景中。

## Ariadne 的特点

- Tree-native：同一个 `tree/` 目录贯穿 discovery、classification 和 phylogeny。
- Feature-space 解释能力：候选序列会被投影到 TPS HMM score space，用于最近参考分配和可视化筛查。
- 面向 CeeSs 的优先识别：如果提供 `TPS/TPS.xlsx` 并安装 ESM 依赖，Ariadne 可以继续筛选更像 `cembrene A / cembrene B` 的 coral-like candidates。
- 建树衔接直接：过滤后的候选可以直接进入 MAFFT alignment 和 IQ-TREE 推断。

## 四阶段流程

1. `discovery`  
   用 query HMM 在蛋白输入或转录组推断的 ORF 中搜索 TPS 候选。
2. `filtering`  
   做 coverage 过滤、最短长度过滤、近重复去冗余，并可额外移除已经存在于参考集中的候选。
3. `classification`  
   将候选和参考一起打到 TPS HMM library 上，构建特征矩阵、做 embedding，并给出最近参考来源。
4. `phylogeny`  
   将过滤后的候选和参考序列合并，运行 MAFFT，再用 IQ-TREE 生成最终系统发育树。

如果没有显式传入 `--query-hmm` 或 `--tps-hmm-dir`，Ariadne 会在 `ariadne/hmm/` 有内置资源时优先复用它们；否则再从 `tree/` 自动构建。

## 结构框架图

<p align="center">
  <img src="./docs/assets/overview_pipeline.svg" alt="Ariadne 结构框架图" width="100%">
</p>

<p align="center">
  这张结构图概括了当前版本的软件框架：同一个 tree-native 参考骨架贯穿 discovery、filtering、特征空间 classification、可选的 CeeSs 打分以及下游 phylogeny。
</p>

## 结果预览

<p align="center">
  <img src="./docs/assets/latest_embedding.svg" alt="Ariadne classification embedding 结果预览" width="100%">
</p>

<p align="center">
  当前仓库自带的本地结果快照来自 <code>result/03_classification/</code>。这里展示的是该目录下的 <code>embedding.svg</code>；在这次本地运行中，共有 <code>36</code> 条 coral-like candidate 被打分，其中 <code>7</code> 条被保留为最终 CeeSs 候选。
</p>

## 安装

推荐环境：Python `3.11`，优先使用仓库自带的 Conda 环境。

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

最小运行依赖：

- Python `>= 3.9`
- `mafft`
- `iqtree` 或 `iqtree2`
- `numpy >= 1.24`
- `openpyxl >= 3.1`
- `pyhmmer >= 0.12.0`
- `pyrodigal >= 3.7.0`
- `scikit-learn >= 1.4`

如果需要 CeeSs 的 ESM 打分，再安装：

```bash
pip install -e '.[esm]'
```

## 快速开始

### 蛋白输入的一键完整流程

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

预期的顶层输出结构：

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
├── 04_phylogeny/
└── pipeline_summary.tsv
```

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

### 只做建树

```bash
ariadne phylogeny \
  --candidates results/02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir results_phylogeny/
```

## 可选的 CeeSs 模式

如果仓库中存在 `TPS/TPS.xlsx`，并且已经安装可选的 ESM 依赖，那么 `run` 和 `classify` 在第三阶段会额外执行一层 CeeSs 打分。当前这层判别逻辑位于 `ariadne.ceess`，使用冻结的 ESM embedding 加一个可训练的小头，默认是 MLP，同时保留旧的 logistic regression 以便复现实验。

这层 CeeSs 头现在是一层直接的多分类：

1. 多分类 TPS 头会直接给每条 `coral-like` candidate 预测细粒度 `Type`。
2. `P(CeeSs)` 由 workbook 中被标记为 CeeSs 正类的所有标签概率求和得到；如果 `TPS.xlsx` 里没有显式的 `CeeSs_group`/`CeeSs` 列，Ariadne 才会回退到只把 `cembrene A` 和 `cembrene B` 当作正类。

常见的额外输出包括：

- `ceess_predictions.tsv`（包含 `esm_type_prediction`、`esm_ceess_probability`，以及按类型展开的 `esm_type_probability_*` 字段）
- `ceess_candidates.tsv`
- `ceess_dark_matter.tsv`
- `ceess_candidates.fasta`
- `ceess_embedding.svg`

## 仓库结构

```text
Ariadne/
├── ariadne/                # 核心包与内置 HMM 资源
├── docs/                   # 项目文档
├── fig/                    # logo 与图像资源
├── input/                  # 示例蛋白输入
├── output/                 # 历史示例输出
├── TPS/                    # 可选的珊瑚 TPS 标注表
├── tree/                   # 默认参考 FASTA 集合
├── environment.yml
├── install.sh
└── pyproject.toml
```

## 文档入口

- [Getting Started](./docs/getting-started.md)
- [Method](./docs/method.md)
- [Tutorials](./docs/tutorials.md)
- [CLI Reference](./docs/cli-reference.md)
- [Outputs](./docs/outputs.md)
- [FAQ](./docs/faq.md)

## 引用

```bibtex
@software{jiang2026ariadne,
  author       = {Jiang, Zhaorui},
  title        = {Ariadne: A coral-centered terpene synthase discovery and CeeSs prioritization platform},
  year         = {2026},
  url          = {https://github.com/zhaoruijiang26/Ariadne},
  version      = {1.0.0}
}
```
