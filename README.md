# Ariadne

<p align="center">
  <img src="fig/logo.png" alt="Ariadne logo" width="280"/>
</p>

<p align="center">
  <strong>A protein-first pipeline for TPS discovery, motif-aware clade assignment, and cembrene benchmarking</strong>
</p>

<p align="center">
  <code>discover -> filter -> motif -> classify -> benchmark</code>
</p>

<p align="center">
  <a href="https://www.python.org/">
    <img src="https://img.shields.io/badge/Python-3.9%2B-3776AB?logo=python&logoColor=white" alt="Python >= 3.9"/>
  </a>
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/License-MIT-16a34a.svg" alt="License MIT"/>
  </a>
  <img src="https://img.shields.io/badge/Version-0.1.0-0f766e.svg" alt="Version 0.1.0"/>
  <a href="https://doi.org/10.1038/s41467-023-44497-0">
    <img src="https://img.shields.io/badge/DOI-10.1038%2Fs41467--023--44497--0-0ea5e9.svg" alt="DOI"/>
  </a>
</p>

## 1. Software Overview

Ariadne 面向 TPS（terpene synthase）候选发现与解释，核心思路是：

- 用查询 HMM 从蛋白 FASTA 中抓出候选。
- 用两阶段 motif 逻辑先判断是否为 TPS，再判断是否更像 cembrene synthase。
- 用多类群参考库把候选投影到低维空间，输出带有 `TPS+ / TPS- / cembrene-like` 标记的 3D cluster context。
- 用 benchmark 模块把最终 FASTA 和理想 FASTA 做精确比对与近似比对。

这一步的 3D 可视化参考了 [AFPK_finder](https://github.com/linzhenjian/AFPK_finder/tree/v1.0) 的高层思路：先构建 profile-score feature space，再做低维投影，观察目标序列是否落入目标类群的聚类区域。

![Ariadne framework](fig/ariadne_framework_nature.svg)

### Key Modules

| Module              | Role                                     | Main outputs                                                                                         |
| ------------------- | ---------------------------------------- | ---------------------------------------------------------------------------------------------------- |
| `discovery.py`      | HMM 命中候选蛋白                         | `01_discovery/candidates.protein.faa`, `candidates.hits.tsv`                                         |
| `filtering.py`      | 覆盖度、长度、去冗余过滤                 | `02_filtering/candidates.filtered.faa`, `filter_report.tsv`                                          |
| `motif.py`          | `DDXXD/E` TPS gate + cembrene motif gate | `04_motif/motif_summary.tsv`, `cembrene_candidates.tsv`, `cembrene_candidates.fasta`                 |
| `classification.py` | 参考类群判定与 3D cluster context        | `03_classification/classification.tsv`, `embedding_3d_sections.svg`, `candidate_cluster_context.tsv` |
| `benchmark.py`      | 预测 FASTA vs 理想 FASTA 对比            | `05_benchmark/benchmark_summary.tsv`, `benchmark_comparison.tsv`                                     |

## 2. Installation

### Quick install

```bash
bash install.sh
ariadne --help
```

### Manual install

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

## 3. Usage

### 3.1 One-command benchmark on this repository

下面这组命令就是当前仓库 `input/` 和 `output/` 的推荐复现实验流程。

```bash
ariadne prepare-references \
  --output-dir benchmark_refs_20260319
```

```bash
ariadne run \
  --protein-folder input \
  --seed-alignment Alignment.fasta \
  --reference-dir benchmark_refs_20260319 \
  --reference-alignment "coralTPS (modified)-cembrene.fasta" \
  --allow-center-fallback \
  --expected-fasta output/output.fasta \
  --output-dir benchmark_run_current_v2b
```

如果你只想单独比较两个 FASTA：

```bash
ariadne compare-fasta \
  --predicted-fasta benchmark_run_current_v2b/04_motif/cembrene_candidates.fasta \
  --expected-fasta output/output.fasta \
  --output-dir benchmark_run_current_v2b/manual_compare
```

### 3.2 Strict vs lenient motif mode

默认是严格模式，不启用中心回退：

```bash
ariadne run \
  --protein-folder input \
  --seed-alignment Alignment.fasta \
  --reference-dir benchmark_refs_20260319 \
  --reference-alignment "coralTPS (modified)-cembrene.fasta" \
  --output-dir benchmark_run_strict
```

宽松模式会打开 `--allow-center-fallback`：

```bash
ariadne run \
  --protein-folder input \
  --seed-alignment Alignment.fasta \
  --reference-dir benchmark_refs_20260319 \
  --reference-alignment "coralTPS (modified)-cembrene.fasta" \
  --allow-center-fallback \
  --output-dir benchmark_run_lenient
```

### 3.3 How to read the outputs

最重要的结果文件：

| Path                                              | Meaning                                                                      |
| ------------------------------------------------- | ---------------------------------------------------------------------------- |
| `03_classification/embedding_3d_sections.svg`     | 3D PCA 切面图。背景是参考类群，前景高亮 `TPS+ / TPS- / cembrene-like` 候选。 |
| `03_classification/candidate_cluster_context.tsv` | 每个候选的 `pc1/pc2/pc3`、候选分组、最近参考类群。                           |
| `03_classification/classification.tsv`            | 每个候选的预测来源、最近邻、置信度。                                         |
| `04_motif/motif_summary.tsv`                      | `DDXXD/E` 是否命中、anchor 是否命中、cembrene-like 是否成立。                |
| `04_motif/cembrene_candidates.fasta`              | 最终导出的 cembrene-like FASTA，可直接与理想结果比较。                       |
| `05_benchmark/benchmark_summary.tsv`              | 总体 benchmark 计数。                                                        |
| `05_benchmark/benchmark_comparison.tsv`           | 逐条序列的 exact / best match 结果。                                         |

### 3.4 Current benchmark result in this repository

基于当前代码，在 `input/all_cn_tp.fasta` 上运行 `--allow-center-fallback` 得到：

| Metric                    | Value |
| ------------------------- | ----- |
| discovery hits            | `100` |
| filtered candidates       | `36`  |
| final cembrene candidates | `7`   |
| ideal sequences           | `13`  |
| exact sequence recovery   | `2`   |
| predicted-only sequences  | `5`   |
| expected-only sequences   | `11`  |

更细的 benchmark 结果来自 `benchmark_run_current_v2b/05_benchmark/benchmark_comparison.tsv`：

- exact recovered: `S_TmTC-1-cembreneA`, `TmTC-1-Cembrene_A`
- near-exact recovered: `S_AbTC-2` 的最佳匹配 identity 为 `0.9925`
- 当前候选整体更偏向 `coral`，但仍有一部分落在 `insect` 近邻区域

当前 3D 聚类上下文统计：

- `candidate:cembrene_like -> coral`: `6`
- `candidate:cembrene_like -> insect`: `1`
- `candidate:tps_positive -> coral`: `20`
- `candidate:tps_positive -> insect`: `5`
- `candidate:tps_negative -> insect`: `4`

这说明现在的流程已经能够：

- 先用 `DDXXD/E` 把 `TPS+` 和 `TPS-` 分开。
- 在 3D 投影里把 cembrene-like 候选作为单独前景标记。
- 把最终导出的 FASTA 和理想结果做可追踪的序列级 benchmark。

### 3.5 Notes

- 当前 benchmark 使用的是仓库内置 `ariadne/tps_hmm`，CLI 会提示这是 legacy library。论文级分析建议传入你自己的 `--tps-hmm-dir`。
- `output/output.fasta` 里存在带 gap 的对齐序列，benchmark 模块会先 ungap 再比较。
- `classification` 现在在 `run` 中会读取 motif 结果，因此 3D 图不再只是“候选 vs 参考”，而是“`TPS+ / TPS- / cembrene-like` vs 多类群参考”。

## 4. Citation

### Software

```text
Ariadne TPS Pipeline (v0.1.0). Accessed YYYY-MM-DD.
```

### Related references

- AFPK_finder repository: [https://github.com/linzhenjian/AFPK_finder/tree/v1.0](https://github.com/linzhenjian/AFPK_finder/tree/v1.0)
- Nature Communications paper: [https://www.nature.com/articles/s41467-023-44497-0](https://www.nature.com/articles/s41467-023-44497-0)
- DOI: `10.1038/s41467-023-44497-0`
