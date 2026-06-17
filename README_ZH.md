<p align="center">
  <img src="docs/assets/ariadne_icon.svg" alt="Ariadne" width="120">
</p>

<h1 align="center">Ariadne</h1>

<p align="center">
  <strong>面向珊瑚萜烯合成酶发现与 cembrene 类合成酶（CeeSs）优先识别的 tree-native 平台</strong>
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
  <img alt="ESM2" src="https://img.shields.io/badge/Optional-ESM2%20CeeSs-2563eb?style=flat-square">
  <img alt="MAFFT + IQ-TREE" src="https://img.shields.io/badge/Phylogeny-MAFFT%20%2B%20IQ--TREE-334155?style=flat-square">
  <a href="./docs/index.md"><img alt="Docs" src="https://img.shields.io/badge/Docs-Read%20the%20Docs-1d4ed8?style=flat-square&logo=readthedocs&logoColor=white"></a>
</p>

---

## 摘要

萜烯合成酶（TPS）产生了自然界中最庞大、结构最多样的天然产物家族。通过基因组挖掘发现候选 TPS 基因已经相对成熟，真正困难的是**产物归属**——判断众多候选合成酶中究竟哪一个负责合成某个特定代谢产物。这个问题在珊瑚（刺胞动物）中尤为突出：它们的基因组编码了数百个候选 TPS，其中作为化学防御核心的 cembrene 类二萜化合物却鲜有经过实验验证的合成酶。

**Ariadne** 是一个 tree-native、四阶段的平台，用于珊瑚 TPS 挖掘以及 **cembrene 类合成酶（CeeSs）** 的优先识别。它从转录组或预测蛋白质组出发，依次执行 profile-HMM 引导的发现、覆盖度与长度过滤、在 TPS profile-HMM 特征空间中的有监督降维分类，以及一个可选的、基于 ESM2 的打分层，用 `P(CeeSs)` 对 coral-like 候选进行排序。单一且经过整理的参考目录 `tree/` 贯穿所有阶段，使候选发现、特征空间解释与最终的最大似然系统发育树始终处于同一个一致的生物学背景中。

## 主要贡献

- **一个参考骨架，四个阶段。** 单一的 `tree/` 目录驱动 discovery、classification 与 phylogeny，消除了各阶段使用不同参考集时产生的人工对接与生物学漂移。
- **可解释的特征空间。** 候选序列在多类群 TPS profile-HMM 打分空间中被嵌入（有监督 LDA，PCA 作为回退），并通过 *k* 近邻投票完成最近参考的标签转移——给出的是几何、邻居与标签，而非一个不透明的分数。
- **基于蛋白质语言模型的 CeeSs 打分。** 当提供标注数据（`TPS/TPS.xlsx`）并安装可选 ESM 依赖时，冻结的 ESM2 主干配合一个轻量可训练头（MLP、逻辑回归，或 Barlow Twins 对比学习变体）对每个 coral-like 候选进行 cembrene A / B 打分。
- **可复现、可审计的输出。** 每次运行都会产出 SVG 嵌入图与系统发育树图、MAFFT 比对、IQ-TREE 系统发育树，以及一条从原始输入到最终候选清单的完整 TSV 审计链。

## 方法

<p align="center">
  <img src="docs/assets/overview_pipeline.svg" alt="Ariadne 流程总览" width="100%">
</p>

Ariadne 是一个顺序执行的四阶段流程，每个阶段也都可作为独立命令调用。

| 阶段 | 命令 | 输入 | 关键输出 |
|---|---|---|---|
| 1 · 发现 | `ariadne discover` | 转录组 / 蛋白 FASTA | `candidates.protein.faa` |
| 2 · 过滤 | `ariadne filter` | 候选 FASTA | `candidates.filtered.faa` |
| 3 · 分类 | `ariadne classify` | 过滤后 FASTA + `tree/` | `classification.tsv`、`embedding.svg` |
| 4 · 建树 | `ariadne phylogeny` | 过滤后 FASTA + `tree/` | `iqtree.treefile`、`phylogeny_preview.svg` |

1. **发现。** 用 Pyrodigal（meta 模式）从转录组装配中预测 ORF，并用基于珊瑚参考比对构建的 profile HMM 搜索翻译后的蛋白；直接提供蛋白 FASTA 时跳过 ORF 预测。
2. **过滤。** 按覆盖度（默认 ≥ 10×）与最短长度（默认 ≥ 300 aa）过滤，并用有界编辑距离在 95% 同一性下折叠近重复序列；与参考匹配的候选会**保留**在 `candidates.filtered.faa` 中，并记录在 `reference_matches.tsv`，从而不会误删已知珊瑚 TPS 基因的新等位变体。
3. **分类。** 所有参考与候选都对 TPS HMM 库打分形成特征向量；矩阵经 z-score 标准化后通过有监督 LDA（对庞大的珊瑚参考集做 *k*-means 子聚类）投影到 3D，PCA 作为回退；再用 *k* 近邻投票给出最近参考标签。当 `TPS/TPS.xlsx` 与 ESM 依赖就绪时，对 coral-like 候选执行 ESM2 子阶段并报告 `P(CeeSs)`。
4. **建树。** 将过滤后的候选与参考合并去冗余，用 MAFFT 比对，再用 IQ-TREE 推断最大似然树，并从生成的 Newick 直接渲染一张紧凑的 SVG 预览图。

## 结果

<p align="center">
  <img src="docs/assets/latest_embedding.svg" alt="TPS 特征空间嵌入结果" width="100%">
</p>

<p align="center">
  <em>某次代表性运行的 TPS profile-HMM 特征空间嵌入：发现 100 个候选 → 过滤后保留 36 个 → 36 个被分类为 coral-like → 在 <code>P(CeeSs) ≥ 0.9</code> 下入选 5 个 CeeSs 候选（ESM2-650M，MLP 头）。</em>
</p>

## 安装

**推荐** —— Python 3.11 + 仓库自带的 Conda 环境：

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

**最小安装** —— 仅核心依赖：

```bash
pip install -e .
# 建树阶段另需 PATH 中存在 mafft 和 iqtree（或 iqtree2）
```

**启用 ESM2 CeeSs 打分：**

```bash
pip install -e '.[esm]'   # 额外安装 torch、transformers、tqdm
```

核心依赖：`numpy ≥ 1.24`、`pyhmmer ≥ 0.12.0`、`pyrodigal ≥ 3.7.0`、`scikit-learn ≥ 1.4`、`openpyxl ≥ 3.1`。

## 快速开始

### 蛋白输入的一键完整流程

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

输出结构：

```text
results/
├── 01_discovery/          # HMM 命中、各样本蛋白 FASTA
├── 02_filtering/          # 过滤后 FASTA、filter_report.tsv、dedupe_clusters.tsv
├── 03_classification/     # classification.tsv、embedding.svg、各候选局部树
│                          #   （存在 TPS/TPS.xlsx 时附带 ceess_* 输出）
├── 04_phylogeny/          # iqtree.treefile、phylogeny_preview.svg
└── pipeline_summary.tsv
```

### 转录组模式

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results_from_transcriptomes/
```

### 单阶段运行

```bash
ariadne classify  --candidates results/02_filtering/candidates.filtered.faa --reference-dir tree/ --output-dir results/03_classification/
ariadne phylogeny --candidates results/02_filtering/candidates.filtered.faa --reference-dir tree/ --output-dir results/04_phylogeny/
```

### 从 Python / Notebook 调用

流程各阶段已在包顶层重新导出，可直接从脚本或 notebook 驱动；裸 `import ariadne` 保持轻量，重型可选依赖在首次使用时才惰性导入（PEP 562）。

```python
import ariadne as ad

discovery = ad.discover_candidates_from_proteins(
    ad.collect_protein_files("input/"), "ariadne/hmm/query.hmm", "results/01_discovery",
)
filtering = ad.filter_candidates(
    discovery["candidate_proteins"], "results/02_filtering", reference_dir="tree/",
)
classification = ad.classify_candidates(
    filtering["filtered_fasta"], reference_dir="tree/",
    output_dir="results/03_classification",
    hmm_dir="ariadne/hmm",   # 内置 TPS HMM 库
    ceess_xlsx=None,         # 设为 "TPS/TPS.xlsx" 以启用 ESM2 CeeSs 打分
)
```

完整可运行示例见 [`examples/tutorial.ipynb`](examples/tutorial.ipynb)。

## CeeSs 打分（可选）

当仓库中存在 `TPS/TPS.xlsx` 且已安装可选 ESM 依赖时，`run` 与 `classify` 会在 HMM 分类后追加一层 ESM2 打分：

1. 从 `TPS.xlsx` 读取标注序列（`Name` / `Protein` / `Type` / `Species`）；
2. 对训练序列与 coral-like 候选计算冻结、均值池化的 ESM2 embedding；
3. 在训练 embedding 上训练一个轻量分类头；
4. 对每个候选打分，`P(CeeSs)` 为所有 CeeSs 正类标签概率之和；
5. 高于 `--ceess-threshold`（默认 0.9）的候选写入 `ceess_candidates.{tsv,fasta}`。

分类头可选 `mlp`（默认）、`logreg` 或 `contrastive`（Barlow Twins）。完整输出说明见 [docs/esm-type.md](./docs/esm-type.md)。

## 命令行

| 命令 | 用途 |
|---|---|
| `ariadne run` | 完整四阶段流程 |
| `ariadne discover` | 第一阶段：HMM 候选发现 |
| `ariadne filter` | 第二阶段：质量过滤与去冗余 |
| `ariadne classify` | 第三阶段：特征空间分类 |
| `ariadne phylogeny` | 第四阶段：MAFFT 比对 + IQ-TREE 建树 |
| `ariadne prepare-references` | 从源数据准备参考 FASTA |
| `ariadne build-hmm` | 从比对构建单个 profile HMM |
| `ariadne build-tps-hmm-library` | 从多个比对构建 TPS HMM 库 |

> 全局参数（`--verbose`、`--log-file`）必须放在子命令之前：`ariadne --verbose run …`

每个命令的完整参数表见 **[CLI Reference](./docs/cli-reference.md)**。

## 仓库结构

```text
Ariadne/
├── ariadne/           # 核心包（含内置 HMM 资源）
├── docs/              # 文档站点（MkDocs + Material）
├── examples/          # 可运行示例脚本与教程 notebook
├── input/             # 示例蛋白输入
├── tree/              # 默认多类群参考 FASTA 集合
├── TPS/               # 用于 CeeSs 打分的珊瑚 TPS 标注表（TPS.xlsx）
├── environment.yml
└── pyproject.toml
```

## 文档

完整文档站点（MkDocs + Material）位于 [`docs/`](./docs/index.md)：

- [Getting Started](./docs/getting-started.md) —— 安装与首次运行
- [Method](./docs/method.md) —— 四阶段设计逐步解析
- [Tutorials](./docs/tutorials.md) —— 实用分析路径
- [CLI Reference](./docs/cli-reference.md) —— 全部命令与参数
- [Outputs](./docs/outputs.md) —— 流程产出的全部文件
- [CeeSs Classifier](./docs/esm-type.md) —— ESM2 打分层
- [Citation](./docs/citation.md)

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
