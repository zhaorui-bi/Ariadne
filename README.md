# Ariadne

<p align="center">
  <img src="fig/ariadne_cvpr_logo.svg" alt="Ariadne CVPR-style logo" width="980">
</p>

<!-- <p align="center">
  <img src="fig/ariadne_icon_minimal.svg" alt="Ariadne minimal icon" width="92">
</p> -->
<p align="center">
  <a href="./README_EN.md">English</a>
</p>

<p align="center">
  🧬 <strong>珊瑚 TPS / CeSS 定向挖掘平台</strong><br>
  🔬 HMM 引导挖掘 · 🧠 motif 校准 · 🌊 珊瑚感知 embedding · 🌳 系统发育分析
</p>

<p align="center">
  <img alt="Python" src="https://img.shields.io/badge/Python-3.11%2B-0F172A?style=flat-square&logo=python&logoColor=white">
  <img alt="MAFFT" src="https://img.shields.io/badge/MAFFT-required-0F766E?style=flat-square">
  <img alt="IQ-TREE" src="https://img.shields.io/badge/IQ--TREE-required-1D4ED8?style=flat-square">
  <img alt="CeSS" src="https://img.shields.io/badge/Target-CeSS%20Mining-F59E0B?style=flat-square">
  <img alt="Benchmark" src="https://img.shields.io/badge/Benchmark-Opt--in-7C3AED?style=flat-square">
</p>

> 🧵 Ariadne 的目标不是只找到“一般 TPS candidate”，而是把 coral 相关数据中的 **CeSS-like / cembrene synthase** 候选收敛成一个更适合实验验证的小集合。

## ✨ 软件功能概览

Ariadne 是一个面向 coral-centered terpene synthase discovery 的命令行平台，当前流程包括：

1. `discovery`
   用 `tree/` 参考自动建 HMM，并从蛋白或 transcriptome 中发现候选 TPS
2. `filtering`
   做 coverage、长度和去冗余过滤
3. `motif`
   做 TPS motif gate，并结合 cembrene family 与 validated CeSS 参考判断 `CeSS-like`
4. `classification`
   在 TPS HMM feature space 中做最近邻分类、clade-aware embedding 和 context tree
5. `phylogeny`
   用 MAFFT + IQ-TREE 构建系统发育树
6. `benchmark`
   与 expected FASTA 对比

## 🚀 当前默认行为

当前版本默认行为已经调整成更适合日常使用：

- ✅ 默认 **不打开 benchmark**
- ✅ 默认 **开启 center fallback**
- ✅ 默认从 `tree/` 自动构建：
  - `query.hmm`
  - TPS HMM library
  - motif 参考
  - phylogeny 参考
- ✅ 默认保留 CeSS 定向挖掘的 calibrated motif 参数
- ✅ benchmark 改成 **显式 opt-in**

这意味着一般情况下你只需要跑：

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

而不是一上来就进入 benchmark 模式。

## 🗂️ 仓库结构

```text
Ariadne/
├── ariadne/                # 核心代码
├── input/                  # 待分析输入
├── output/                 # 预期输出 / benchmark 示例
├── tree/                   # 多物种 TPS 参考，同时也是默认 HMM 来源
├── fig/                    # logo / figures
├── environment.yml         # conda 环境
├── install.sh
└── pyproject.toml
```

### 📁 数据约定

- `input/`
  日常分析输入目录。当前仓库示例是 `input/all_tps.fasta`
- `tree/`
  当前最关键的参考目录，既用于 HMM 构建，也用于 motif / classification / phylogeny
- `output/`
  预期输出示例目录。当前仓库保留文件名 `output/outpu.fasta`

## 🧩 核心模块

| 模块                        | 作用                                                                                |
| --------------------------- | ----------------------------------------------------------------------------------- |
| `ariadne/cli.py`            | CLI 入口，定义 `run`、`discover`、`motif`、`classify`、`phylogeny`、`compare-fasta` |
| `ariadne/discovery.py`      | HMM 构建、候选发现、ORF 预测                                                        |
| `ariadne/filtering.py`      | coverage、长度、去冗余过滤                                                          |
| `ariadne/motif.py`          | TPS motif gate、CeSS 判断、priority ranking                                         |
| `ariadne/classification.py` | feature matrix、embedding、最近邻分类、context tree                                 |
| `ariadne/phylogeny.py`      | MAFFT + IQ-TREE                                                                     |
| `ariadne/benchmark.py`      | FASTA benchmark 对比                                                                |
| `ariadne/fasta_utils.py`    | FASTA / TSV 读写与清洗                                                              |

## 🛠️ 安装

### 1. 推荐：conda

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
```

### 2. venv

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
bash install.sh
```

### 3. 手动安装

```bash
git clone https://github.com/zhaoruijiang26/Ariadne.git
cd Ariadne
python -m venv .venv
source .venv/bin/activate
python -m pip install -U pip
python -m pip install -e .
```

### 📦 运行依赖

- Python `>= 3.9`
- `mafft`
- `iqtree`
- `numpy >= 1.24`
- `openpyxl >= 3.1`
- `pyhmmer >= 0.12.0`
- `pyrodigal >= 3.7.0`
- `scikit-learn >= 1.4`

推荐 Python `3.11`。

## ⚙️ 默认参数画像

下面这些是当前更适合日常 coral / CeSS mining 的默认设置：

| 参数                                | 当前默认值 | 含义                                                    |
| ----------------------------------- | ---------- | ------------------------------------------------------- |
| Benchmark                           | `off`      | 只有显式传 `--enable-benchmark` 才生成 `05_benchmark/`  |
| Center fallback                     | `on`       | 缺失 direct anchor 时，默认回退到中心窗口继续 CeSS 判断 |
| `validated_cess_identity_threshold` | `0.95`     | 高置信 validated CeSS 支持阈值                          |
| `cembrene_identity_threshold`       | `0.75`     | cembrene family 级支持阈值                              |
| `cembrene_margin_threshold`         | `0.10`     | 相对 non-cembrene 的最小 identity margin                |
| IQ-TREE fast mode                   | `on`       | 默认使用 `LG + --fast` 跑完整流程                       |

## 🧪 快速开始

### 1. 日常挖掘模式

这是最推荐的日常入口。✅ 不 benchmark，✅ 默认 center fallback，✅ 自动从 `tree/` 建模：

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --output-dir results/
```

### 2. CeSS 校准模式

如果你手头有已验证 CeSS FASTA，建议打开校准：

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --reference-alignment tree/coral.fasta \
  --validated-cess-fasta output/outpu.fasta \
  --output-dir results_calibrated/
```

### 3. Benchmark 模式

只有需要做方法对照时再打开 benchmark：

```bash
ariadne run \
  --protein-folder input/ \
  --reference-dir tree/ \
  --reference-alignment tree/coral.fasta \
  --validated-cess-fasta output/outpu.fasta \
  --expected-fasta output/outpu.fasta \
  --enable-benchmark \
  --output-dir results_benchmark/
```

### 4. 转录组输入

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir tree/ \
  --output-dir results_from_transcriptome/
```

## 🔍 分阶段使用

### `discover`

```bash
ariadne discover \
  --protein-folder input/ \
  --hmm query.hmm \
  --output-dir 01_discovery/
```

### `filter`

```bash
ariadne filter \
  --input-fasta 01_discovery/candidates.protein.faa \
  --output-dir 02_filtering/
```

### `motif`

```bash
ariadne motif \
  --candidates 02_filtering/candidates.filtered.faa \
  --reference-alignment tree/coral.fasta \
  --validated-cess-fasta output/outpu.fasta \
  --output-dir 04_motif/
```

常用可调参数：

- `--disable-center-fallback`
  关闭默认的 center fallback，切回更严格模式
- `--validated-cess-identity-threshold`
  调整 validated CeSS 支持强度
- `--cembrene-identity-threshold`
  调整 family-level identity 下限
- `--cembrene-margin-threshold`
  调整与 non-cembrene 的分离 margin

### `classify`

```bash
ariadne classify \
  --candidates 02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir 03_classification/
```

### `phylogeny`

```bash
ariadne phylogeny \
  --candidates 02_filtering/candidates.filtered.faa \
  --reference-dir tree/ \
  --output-dir 06_phylogeny/
```

### `compare-fasta`

```bash
ariadne compare-fasta \
  --predicted-fasta 04_motif/cess_candidates.fasta \
  --expected-fasta output/outpu.fasta \
  --output-dir 05_benchmark/
```

## 📤 输出目录地图

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
│   ├── classification.tsv
│   ├── cess_priority_ranking.tsv
│   ├── embedding.tsv
│   ├── embedding.svg
│   ├── embedding_3d_sections.svg
│   ├── global_context_tree.nwk
│   └── _auto_tps_hmms/
├── 04_motif/
│   ├── motif_summary.tsv
│   ├── targeted_mining_summary.tsv
│   ├── cess_candidates.tsv
│   ├── cess_candidates.fasta
│   ├── cess_priority_ranking.tsv
│   └── motif_windows.svg
├── 05_benchmark/          # 只有 --enable-benchmark 才会生成
└── 06_phylogeny/
```

### 🎯 最值得先看的结果

- `04_motif/targeted_mining_summary.tsv`
  看整体候选数量、TPS+/TPS-、CeSS-like 数量
- `04_motif/cess_candidates.fasta`
  最终 CeSS-like 候选序列
- `03_classification/cess_priority_ranking.tsv`
  最适合做实验验证优先级筛选的表
- `03_classification/embedding.svg`
  当前的 clade-aware embedding 总览
- `06_phylogeny/iqtree.treefile`
  系统发育树结果

## 🧠 Embedding 设计

当前 `classification` 不是单纯画一个“拥挤的 PCA 点图”，而是优先使用：

- `lda_coral_subclade_spread`
  先按 reference clades / candidate groups 做监督式分离
- 对 `coral` 再进一步拆成：
  - `cembrene_a`
  - `cembrene_b`
  - `cembrene_c`
  - 多个 `tps_subclade_*`

所以新版本的 `embedding.svg` 会更适合直接看：

- reference clades 是否分开
- coral 内部 CeSS/TPS 子簇是否分开
- candidate TPS+/TPS-/CeSS-like 分布在哪里

## 🧵 CeSS 判定逻辑

当前 CeSS 判定逻辑如下：

- 如果与 validated CeSS motif window identity `>= 0.95`
  直接支持 `CeSS-like`
- 如果 direct anchor 命中，同时：
  - 对 cembrene 参考 identity `>= 0.75`
  - 且相对 non-cembrene 的 margin `>= 0.10`
    也支持 `CeSS-like`
- 如果只是 center fallback 命中但缺乏 validated CeSS 支持
  不直接判成 CeSS-like，但会进入 `cess_priority_ranking.tsv`

## 📊 仓库内实跑结果

当前仓库已经实跑出一套更完整的默认结果：

- 📂 [tmp_run_calibrated_v3/](tmp_run_calibrated_v3)
- 🖼️ [tmp_run_calibrated_v3/03_classification/embedding.svg](tmp_run_calibrated_v3/03_classification/embedding.svg)
- 🏁 [tmp_run_calibrated_v3/05_benchmark/benchmark_summary.tsv](tmp_run_calibrated_v3/05_benchmark/benchmark_summary.tsv)

关键摘要：

- `36` input candidates
- `32` TPS-positive
- `4` TPS-negative
- `4` CeSS-like
- benchmark 模式下 `3` exact matches
- 另有 `1` 条高相似近似命中

## 🌳 HMM 与系统树说明

- `tree/` 现在是默认参考入口，不再以 `Alignment.fasta` 为主入口
- 如果 `tree/*.fasta` 不是 MSA，Ariadne 会先自动调用 MAFFT，再做 `hmmbuild`
- `run` 与 `phylogeny` 默认使用 IQ-TREE 的 `LG + --fast`
- 如果你只想挖掘，不想 benchmark，不需要提供 `--expected-fasta`

## 🧰 旧脚本兼容性

这些旧脚本还在仓库中保留，但新的主流程优先使用 CLI：

- `ariadne/filter_contigs.py`
- `ariadne/filter_contigs_long.py`
- `ariadne/filter_coverage.py`
- `ariadne/DupRemover.py`
- `ariadne/visualization.py`

推荐优先使用：

- `ariadne run`
- `ariadne discover`
- `ariadne filter`
- `ariadne motif`
- `ariadne classify`
- `ariadne phylogeny`
- `ariadne compare-fasta`

## 📎 软件定位

如果你把 Ariadne 当成论文中的平台来理解，可以把它看成：

> 🪸 一个以 coral TPS / CeSS 为中心、偏实验筛选友好的 discovery + prioritization + phylogeny workflow。

## 🎨 Logo 资产

- `fig/ariadne_cvpr_logo.svg`
  横版展示 logo，适合 README banner、项目介绍页、幻灯片封面
- `fig/ariadne_icon_minimal.svg`
  极简 icon，适合论文图角标、GitHub avatar、favicon 风格展示

## License

MIT License
