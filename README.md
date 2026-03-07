# Ariadne

Ariadne 现在被整理成一个可直接运行的四模块动物 TPS 挖掘 pipeline，并保留了旧脚本兼容入口。

## 这次重构做了什么

- 用 `python -m ariadne` 提供统一 CLI，覆盖四个模块和端到端运行。
- 用 `pyrodigal` 替代缺失的 `prodigal` 二进制，用 `pyhmmer` 替代 discovery HMM 搜索。
- 用仓库内 `AFLP_finder-main/hmm/*.hmm` 做第 3 模块的 AFLP 风格特征提取，不再依赖 R。
- 用纯 Python 重写过滤、去重和 motif 分析，去掉旧脚本对 Biopython 的隐式依赖。
- 增加 `prepare-demo` 和新的 `code.sh`，可以直接做 smoke test。
- 整理珊瑚 FASTA 与昆虫 Excel 参考数据的预处理流程。

## 当前仓库结构

- `ariadne/cli.py`：主命令行入口。
- `ariadne/discovery.py`：模块 1，转录组 ORF 预测 + HMM 搜索。
- `ariadne/filtering.py`：模块 2，coverage / length / 95% 相似性过滤。
- `ariadne/classification.py`：模块 3，AFLP 风格 HMM 特征分类、局部树和 embedding。
- `ariadne/motif.py`：模块 4，210 aa 附近 motif window 分析与 cembrene 候选判断。
- `ariadne/references.py`：珊瑚 FASTA、昆虫 Excel 和额外参考集整理。
- `ariadne/demo.py`：构建一个可完整跑通的 demo workspace。
- `AFLP_finder-main/`：保留原 AFLP 资料和 30 个 HMM profile。
- `code.sh`：一键 demo / smoke test。

## 安装

推荐用本地虚拟环境：

```bash
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

依赖在 `pyproject.toml` 中声明：

- `pyrodigal`
- `pyhmmer`
- `openpyxl`
- `numpy`

## 四个模块对应关系

### 模块 1：转录组预测与 HMM 搜索

- 输入：转录组 FASTA
- 步骤：`pyrodigal` 预测 ORF → 翻译蛋白 → `pyhmmer` HMM 搜索 → 抽取 hits
- 输出：
  - `all_predicted_proteins.faa`
  - `candidates.protein.faa`
  - `candidates.orf.fna`
  - `candidates.hits.tsv`

### 模块 2：序列过滤与质量控制

- 默认规则：
  - `coverage >= 10`
  - `length >= 300 aa`
  - `identity_threshold = 0.95`
- 输出：
  - `candidates.filtered.faa`
  - `filter_report.tsv`
  - `dedupe_clusters.tsv`
  - `manual_review.tsv`

说明：

- 旧脚本里的 `500 aa` 阈值不是硬编码了，现在通过 `--min-length 500` 即可复现。
- `manual_review.tsv` 会把需要人工检查的序列列出来，方便你继续做视觉检查。

### 模块 3：系统发育 / AFLP 风格结构分类

- 使用 `AFLP_finder-main/hmm/*.hmm` 的 30 个 profile 为每条序列打分。
- 生成 AFLP 风格特征矩阵，不需要 R / t-SNE 也能跑。
- 输出：
  - `aflp_features.tsv`
  - `embedding.tsv`
  - `embedding.svg`
  - `classification.tsv`
  - `nearest_neighbors.tsv`
  - `trees/*.nwk`

这里的分类逻辑是：

- 先把参考序列和候选序列都映射到同一个 HMM feature space；
- 再用最近邻判断候选更靠近哪个来源类群；
- 同时导出局部 UPGMA tree 方便后续看 candidate 周围的邻居。

### 模块 4：功能基序与 cembrene specificity 分析

- 默认先搜 `CFDVL.` anchor；
- 如果参考或候选中找不到完全匹配，就回退到 `210 aa` 附近窗口；
- 再把候选和珊瑚 reference 中标记为 `cembrene` 的序列做窗口相似度比较。

输出：

- `motif_summary.tsv`
- `motif_windows.fasta`
- `motif_windows.svg`
- `cembrene_candidates.tsv`

## 参考库准备

仓库当前直接带了两类参考：

- 珊瑚：`ariadne/coralTPS (modified)-cembrene.fasta`
- 昆虫：`ariadne/Insecta TPS.xlsx`

可以先整理成统一 reference 目录：

```bash
python -m ariadne prepare-references \
  --output-dir references_clean
```

如果你后续补齐 bacteria / fungal / plant FASTA，可以继续追加：

```bash
python -m ariadne prepare-references \
  --output-dir references_clean \
  --extra-fasta bacteria=/path/to/bacteria.fasta \
  --extra-fasta fungal=/path/to/fungal.fasta \
  --extra-fasta plant=/path/to/plant.fasta
```

## Demo / 测试

### 一键跑通

```bash
bash ./code.sh
```

它会：

- 创建本地虚拟环境；
- 安装当前项目；
- 生成 demo 输入；
- 跑完整的 4 模块；
- 输出 `demo_workspace/demo_output/pipeline_summary.tsv`

### 手动 smoke test

```bash
python -m ariadne prepare-demo --output-dir /tmp/ariadne_demo

python -m ariadne run \
  --transcriptomes /tmp/ariadne_demo/transcriptomes/demo_sample_transcripts.fasta \
  --seed-alignment /tmp/ariadne_demo/seed_alignment.fasta \
  --reference-dir /tmp/ariadne_demo/references \
  --output-dir /tmp/ariadne_demo_run
```

这会生成：

- discovery hits：1 条完整 candidate + 1 条低覆盖/短序列被过滤掉
- classification：candidate 被判为 `coral`
- motif：candidate 被判为 `cembrene_like = yes`

## 实际数据运行

假设你已经有：

- 一个对齐好的 seed alignment：`seed_alignment.fasta`
- 一个整理好的参考目录：`references_clean/`
- 多个转录组：`assembly/sample1_transcripts.fasta` 等

可以直接跑：

```bash
python -m ariadne run \
  --transcriptomes assembly/sample1_transcripts.fasta assembly/sample2_transcripts.fasta \
  --seed-alignment seed_alignment.fasta \
  --reference-dir references_clean \
  --output-dir results_run_01
```

常用参数：

- `--min-coverage 10`
- `--min-length 300`
- `--min-length 500`
- `--identity-threshold 0.95`
- `--top-k 5`
- `--tree-neighbors 12`

## 分模块运行

### 1. 构建 query HMM

```bash
python -m ariadne build-hmm \
  --alignment seed_alignment.fasta \
  --output query.hmm
```

### 2. discovery

```bash
python -m ariadne discover \
  --transcriptomes assembly/sample1_transcripts.fasta \
  --hmm query.hmm \
  --output-dir step1_discovery
```

### 3. filtering

```bash
python -m ariadne filter \
  --input-fasta step1_discovery/candidates.protein.faa \
  --output-dir step2_filter \
  --min-coverage 10 \
  --min-length 300 \
  --identity-threshold 0.95
```

### 4. classification

```bash
python -m ariadne classify \
  --candidates step2_filter/candidates.filtered.faa \
  --reference-dir references_clean \
  --output-dir step3_classify
```

### 5. motif / cembrene specificity

```bash
python -m ariadne motif \
  --candidates step2_filter/candidates.filtered.faa \
  --coral-reference references_clean/coral.fasta \
  --output-dir step4_motif
```

## 旧脚本兼容

以下旧入口已保留，但内部已经改成新的纯 Python 实现：

- `ariadne/filter_coverage.py`
- `ariadne/filter_contigs.py`
- `ariadne/filter_contigs_long.py`
- `ariadne/DupRemover.py`

例如：

```bash
python ariadne/filter_coverage.py 10 input.fasta
python ariadne/filter_contigs.py 500 input.fasta
python ariadne/DupRemover.py -i input.fasta
```

## 目前还需要你补的数据

基于当前 repo，我已经把能自动化的部分都接上了，但以下真实分析数据仍然建议你后续补齐：

- bacteria TPS FASTA
- fungal TPS FASTA
- plant TPS FASTA
- 你实际使用的 seed alignment / HMM training 序列

这些数据一旦补齐，只需要通过 `prepare-references` 或 `--extra-fasta` 加进去即可。

## 已完成的本地验证

我已经在当前仓库里完成了 smoke test：

```bash
python -m ariadne prepare-demo --output-dir /tmp/ariadne_demo_smoke4
python -m ariadne run \
  --transcriptomes /tmp/ariadne_demo_smoke4/transcriptomes/demo_sample_transcripts.fasta \
  --seed-alignment /tmp/ariadne_demo_smoke4/seed_alignment.fasta \
  --reference-dir /tmp/ariadne_demo_smoke4/references \
  --output-dir /tmp/ariadne_demo_run4
```

关键验证结果：

- `01_discovery/candidates.hits.tsv`：命中 2 条候选
- `02_filtering/filter_report.tsv`：低 coverage / 短序列被剔除，只保留 1 条完整候选
- `03_classification/classification.tsv`：候选被分类为 `coral`
- `04_motif/motif_summary.tsv`：候选被标记为 `predicted_cembrene_like = yes`

## 备注

- 仓库里的 `AFLP_finder-main/` 没被删除，仍然保留作参考资料和 HMM profile 来源。
- 现在默认流程不依赖外部 `prodigal`、`hmmsearch`、`R`、`mafft`、`iqtree`、`meme`。
- 如果你后续要做发表级系统树，仍然可以拿 `candidates.filtered.faa` 和整理后的 references 再单独跑 `mafft + iqtree`。
