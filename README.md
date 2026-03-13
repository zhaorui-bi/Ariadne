# Ariadne

<p align="center">
  <img src="fig/logo.png" alt="Ariadne Logo" width="320"/>
</p>

<p align="center">
  <strong>Protein-first TPS discovery, clade assignment, and cembrene motif analysis pipeline</strong>
</p>

<p align="center">
  <code>discover -> filter -> classify -> motif</code>
</p>

<p align="center">
  <a href="https://www.python.org/">
    <img src="https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white" alt="Python >= 3.9"/>
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-22c55e.svg" alt="License: MIT"/>
  </a>
  <img src="https://img.shields.io/badge/Version-0.1.0-7c3aed.svg" alt="Version 0.1.0"/>
  <a href="https://doi.org/10.1038/s41467-023-44497-0">
    <img src="https://img.shields.io/badge/DOI-10.1038%2Fs41467--023--44497--0-0ea5e9.svg" alt="DOI"/>
  </a>
</p>

---

## 第一部分：软件介绍（Software Overview）

### 1.1 软件定位

Ariadne 是一个面向 TPS（terpene synthase）序列发现与注释的分析管线，核心目标是：

1. 从蛋白或转录组数据中筛选 TPS 候选序列。
2. 将候选序列投影到多类群参考空间（bacteria/fungal/plant/insect/coral）中进行来源判定。
3. 以两阶段 motif 逻辑评估 cembrene 相关特征（支持严格/宽松两种模式）。

### 1.2 核心能力

- 模块化流程：`discover`、`filter`、`classify`、`motif`。
- 多类群参考库自动构建：支持自动识别本地标准文件名。
- 类群证据输出：`classification.tsv`、`assignment_summary.tsv`、`global_context_tree.nwk`、`embedding_3d_sections.svg`。
- 模体判定两种策略：
  - 严格模式（默认）：TPS 通过后必须命中 anchor 再做 cembrene 判断。
  - 宽松模式（`--allow-center-fallback`）：anchor 缺失时允许中心窗口回退判断。

### 1.3 Nature 风格架构图

![Ariadne Nature-style Architecture](fig/ariadne_framework_nature.svg)

**Figure 1.** 该图与当前代码实现对齐，展示了输入层、四模块主线、strict/lenient motif 开关以及用于论文汇报的关键输出。

### 1.4 关键模块与产物

| 模块 | 主要功能 | 主要产物 |
|---|---|---|
| `discovery.py` | HMM 搜索候选 TPS | `candidates.protein.faa`, `candidates.hits.tsv` |
| `filtering.py` | 覆盖度/长度过滤与去冗余 | `candidates.filtered.faa`, `filter_report.tsv` |
| `classification.py` | 多类群判定、嵌入、树构建 | `classification.tsv`, `assignment_summary.tsv`, `global_context_tree.nwk`, `embedding_3d_sections.svg` |
| `motif.py` | TPS gate + cembrene motif 判定 | `motif_summary.tsv`, `cembrene_candidates.tsv`, `motif_windows.svg` |
| `visualization.py` | t-SNE/聚类可视化 | `tsne-db*.tsv`, `tsne*.svg` |

---

## 第二部分：安装（Installation）

### 2.1 环境要求

- Python `>=3.9`
- 建议使用虚拟环境

### 2.2 快速安装

```bash
bash install.sh
ariadne --help
```

### 2.3 手动安装

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

---

## 第三部分：使用说明（Usage）

### 3.0 快速开始一键复制命令块（Quick Copy）

> GitHub 中每个代码块右上角都有复制按钮，可一键复制。

1) 构建参考库（自动识别）

```bash
ariadne prepare-references \
  --output-dir /private/tmp/ariadne_refs_auto
```

2) 构建 TPS HMM 库

```bash
ariadne build-tps-hmm-library \
  --alignment tps_core=Alignment.fasta \
  --output-dir /private/tmp/ariadne_tps_hmm_real
```

3) 运行完整流程（严格模式）

```bash
ariadne run \
  --protein-folder ./protein \
  --seed-alignment ./Alignment.fasta \
  --reference-dir /private/tmp/ariadne_refs_auto \
  --reference-alignment "./coralTPS (modified)-cembrene.fasta" \
  --tps-hmm-dir /private/tmp/ariadne_tps_hmm_real \
  --output-dir /private/tmp/ariadne_run_strict
```

4) 运行完整流程（宽松模式）

```bash
ariadne run \
  --protein-folder ./protein \
  --seed-alignment ./Alignment.fasta \
  --reference-dir /private/tmp/ariadne_refs_auto \
  --reference-alignment "./coralTPS (modified)-cembrene.fasta" \
  --tps-hmm-dir /private/tmp/ariadne_tps_hmm_real \
  --allow-center-fallback \
  --output-dir /private/tmp/ariadne_run_lenient
```

### 3.1 构建多类群参考库（推荐自动模式）

```bash
ariadne prepare-references \
  --output-dir /private/tmp/ariadne_refs_auto
```

自动识别（存在则加载）：

- `coralTPS (modified)-cembrene.fasta`
- `Insecta TPS.xlsx`
- `bacteria.fasta`（也支持 `.fa/.faa`）
- `fungal.fasta` 或 `fungi.fasta`（也支持 `.fa/.faa`）
- `plant.fasta`（也支持 `.fa/.faa`）

手动指定路径示例：

```bash
ariadne prepare-references \
  --coral "./coralTPS (modified)-cembrene.fasta" \
  --insect-xlsx "./Insecta TPS.xlsx" \
  --bacteria-fasta "./bacteria.fasta" \
  --fungi-fasta "./fungi.fasta" \
  --plant-fasta "./plant.fasta" \
  --output-dir /private/tmp/ariadne_refs_manual
```

### 3.2 构建 TPS HMM 库

```bash
ariadne build-tps-hmm-library \
  --alignment tps_core=Alignment.fasta \
  --output-dir /private/tmp/ariadne_tps_hmm_real
```

### 3.3 运行完整流程（严格/宽松）

严格模式（默认）：

```bash
ariadne run \
  --protein-folder ./protein \
  --seed-alignment ./Alignment.fasta \
  --reference-dir /private/tmp/ariadne_refs_auto \
  --reference-alignment "./coralTPS (modified)-cembrene.fasta" \
  --tps-hmm-dir /private/tmp/ariadne_tps_hmm_real \
  --output-dir /private/tmp/ariadne_run_strict
```

宽松模式（启用中心窗口回退）：

```bash
ariadne run \
  --protein-folder ./protein \
  --seed-alignment ./Alignment.fasta \
  --reference-dir /private/tmp/ariadne_refs_auto \
  --reference-alignment "./coralTPS (modified)-cembrene.fasta" \
  --tps-hmm-dir /private/tmp/ariadne_tps_hmm_real \
  --allow-center-fallback \
  --output-dir /private/tmp/ariadne_run_lenient
```

### 3.4 关键结果怎么看

优先查看：

1. `run_output/03_classification/assignment_summary.tsv`（类群汇总）
2. `run_output/03_classification/classification.tsv`（每条候选来源）
3. `run_output/03_classification/global_context_tree.nwk`（全局树）
4. `run_output/03_classification/embedding_3d_sections.svg`（3D 切面）
5. `run_output/04_motif/motif_summary.tsv`（两阶段门控细节）
6. `run_output/04_motif/cembrene_candidates.tsv`（最终 cembrene-like 候选）

### 3.5 当前仓库最新实测（2026-03-13）

参考库：`/private/tmp/ariadne_refs_auto_20260313`

- total: `774`
- insect: `420`
- coral: `240`
- bacteria: `47`
- plant: `41`
- fungal: `26`

严格模式结果（`/private/tmp/ariadne_multiclade_run_auto_strict_20260313/run_output`）：

- classification: `36` candidates (`coral=36`, mean confidence `1.0`)
- motif: `is_tps=yes 32`, `is_tps=no 4`
- cembrene candidates: `0`

宽松模式结果（`/private/tmp/ariadne_multiclade_run_auto_lenient_20260313/run_output`）：

- motif `predicted_cembrene_like=yes`: `7`
- `cembrene_candidates.tsv`: `8` lines (含表头)

### 3.6 常见排错

- `candidates.hits.tsv` 很少或为空：检查 `--seed-alignment` 和输入 FASTA 质量。
- `filter_report.tsv` 大量 `low_coverage`：调低 `--min-coverage`。
- `filter_report.tsv` 大量 `too_short`：调低 `--min-length`。
- 严格模式下 `cembrene_candidates.tsv` 为空：常见且合理，可切换 `--allow-center-fallback` 对比。

---

## 第四部分：引用（Citation）

### 4.1 软件引用

如果 Ariadne 被用于论文或报告，请至少引用本仓库并注明版本与访问日期：

```text
Ariadne TPS Pipeline (v0.1.0). Available at: <your-repo-url>. Accessed: YYYY-MM-DD.
```

### 4.2 方法学参考文献（你提供的 Nature 链接）

- URL: https://www.nature.com/articles/s41467-023-44497-0
- DOI: `10.1038/s41467-023-44497-0`

BibTeX（可直接放入论文参考文献库）：

```bibtex
@article{s41467_023_44497_0,
  doi = {10.1038/s41467-023-44497-0},
  url = {https://www.nature.com/articles/s41467-023-44497-0},
  journal = {Nature Communications},
  year = {2023}
}
```

> 建议投稿前再次核对该条目的题名、作者列表与研究主题是否与当前稿件完全匹配。
