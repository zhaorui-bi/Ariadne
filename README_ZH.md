<p align="center">
  <img src="docs/logo.png" alt="Ariadne 萜烯合成酶发现平台" width="380">
</p>

<p align="center">
  <strong>珊瑚萜烯合成酶的全基因组定向挖掘平台</strong><br>
  <em>Journal of Natural Products</em>（2026）&nbsp;·&nbsp; DOI: <a href="https://doi.org/10.1021/acs.jnatprod.6c00686">10.1021/acs.jnatprod.6c00686</a>
</p>

<p align="center">
  <a href="./README.md">English</a> &nbsp;·&nbsp;
  <a href="https://ariadne-platform.readthedocs.io/en/latest/">Read the Docs</a> &nbsp;·&nbsp;
  <a href="https://doi.org/10.1021/acs.jnatprod.6c00686">论文</a> &nbsp;·&nbsp;
  <a href="#引用">引用</a>
</p>

<p align="center">
  <a href="https://ariadne-platform.readthedocs.io/en/latest/"><img alt="Read the Docs" src="https://img.shields.io/readthedocs/ariadne-platform/latest?style=flat-square&label=Read%20the%20Docs&color=1d4ed8"></a>
  <a href="https://doi.org/10.1021/acs.jnatprod.6c00686"><img alt="Journal of Natural Products 2026" src="https://img.shields.io/badge/J.%20Nat.%20Prod.-2026-084c61?style=flat-square"></a>
  <a href="https://doi.org/10.1021/acs.jnatprod.6c00686"><img alt="DOI" src="https://img.shields.io/badge/DOI-10.1021%2Facs.jnatprod.6c00686-084c61?style=flat-square"></a>
  <img alt="Python" src="https://img.shields.io/badge/Python-%E2%89%A53.9-0b132b?style=flat-square&logo=python&logoColor=white">
  <img alt="Version" src="https://img.shields.io/badge/Version-1.1.0-0f766e?style=flat-square">
  <img alt="License" src="https://img.shields.io/badge/License-MIT-c2410c?style=flat-square">
</p>

---

**Ariadne** 是 Wang 等在 *J. Nat. Prod.* **2026** 中报道的计算平台，用于珊瑚萜烯合成酶（TPS）的全基因组定向挖掘，并对 **cembrene-class synthase（CeeSs）** 进行优先级排序。它将转录组或预测蛋白 FASTA 转换为可审计的四阶段流程：Discovery、Filtering、Classification 与 PCA/LDA 可视化；在提供 `TPS.xlsx` 并安装 `[esm]` 依赖时，还可对 coral-like 候选进行 ESM2 打分。

> 基因组挖掘是发现天然产物酶的有效策略，但 TPS 搜索常反复得到已知产物。Ariadne 针对珊瑚中 14 元碳环 cembranoid 途径而设计，目标是给出可实验验证的 CeeSs，而不是一份未排序的同源蛋白清单。

## 论文结果

以下为正式发表文章中的结果，不表示本地一次软件运行的输出：

<table>
  <tr>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/flask.svg" width="28" height="28" alt=""><br>
      <strong>5 个 CeeSs</strong><br>
      从珊瑚 TPS 同源蛋白中鉴定的 cembrene 合成酶
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/target.svg" width="28" height="28" alt=""><br>
      <strong>80% 准确率</strong><br>
      90 个 TPS 同源蛋白的实验验证率
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/ring.svg" width="28" height="28" alt=""><br>
      <strong>14 元碳环</strong><br>
      优先挖掘的 cembranoid 骨架
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/genome.svg" width="28" height="28" alt=""><br>
      <strong>多系起源</strong><br>
      进化轨迹用于指导酶工程
    </td>
  </tr>
</table>

发表集合包含两个与已知 cembrene B 合成酶高度相似的酶，以及三个经酵母异源表达验证的低相似度酶。文章中的系统发育分析，以及对早期分歧珊瑚酶的改造（该酶可合成此前未报道的 cembrene 骨架），属于由候选序列推动的后续科学研究，不是本仓库中的额外命令行步骤。

## 科学背景

Cembranoids 是一类以 14 元碳环为骨架的珊瑚二萜。这类骨架的化学合成仍然困难；此前六项独立研究中的十五个珊瑚酶，仅获得两个生物合成前体。常规 TPS 基因组挖掘容易反复发现同一产物家族。Ariadne 将 profile-HMM 发现、质量过滤、参考空间分类和面向 CeeSs 的可视化结合起来，把湿实验精力集中到新颖或表征不足的珊瑚酶上。

## 平台流程

<table>
  <tr>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/discovery.svg" width="28" height="28" alt=""><br>
      <strong>I. Discovery</strong><br>
      对转录组 ORF 或蛋白 FASTA 进行 HMM 搜索
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/filter.svg" width="28" height="28" alt=""><br>
      <strong>II. Filtering</strong><br>
      覆盖度、长度与近重复控制
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/classify.svg" width="28" height="28" alt=""><br>
      <strong>III. Classification</strong><br>
      参考空间分类，并保留近邻证据
    </td>
    <td align="center" valign="top" width="25%">
      <img src="docs/images/icons/visualize.svg" width="28" height="28" alt=""><br>
      <strong>IV. Visualization</strong><br>
      PCA/LDA 投影，用于候选人工筛选
    </td>
  </tr>
</table>

<p align="center">
  <img src="docs/images/algorithm-framework.png" alt="Ariadne 四阶段算法框架" width="880">
</p>

## 核心能力

<table>
  <tr>
    <td valign="top" width="50%">
      <img src="docs/images/icons/genome.svg" width="22" height="22" alt="">
      <strong> 定向珊瑚 TPS 挖掘</strong><br>
      基于用户提供的参考 FASTA 目录进行基因组/转录组筛选；若无预构建 HMM，可由参考序列生成。
    </td>
    <td valign="top" width="50%">
      <img src="docs/images/icons/output.svg" width="22" height="22" alt="">
      <strong> 可审计中间结果</strong><br>
      每个阶段输出 TSV 证据，候选可被保留、剔除或重排，而不是只给出黑箱终点。
    </td>
  </tr>
  <tr>
    <td valign="top">
      <img src="docs/images/icons/target.svg" width="22" height="22" alt="">
      <strong> CeeSs 优先级</strong><br>
      在提供 <code>TPS.xlsx</code> 并安装 <code>[esm]</code> 时，对 coral-like 候选进行可选 ESM2 打分。
    </td>
    <td valign="top">
      <img src="docs/images/icons/visualize.svg" width="22" height="22" alt="">
      <strong> 论文风格图件</strong><br>
      输出三面板 LDA/PCA SVG（LD1–LD2、LD1–LD3、LD2–LD3），供异源表达前检查。
    </td>
  </tr>
</table>

## 安装

推荐方式：

```bash
git clone https://github.com/zhaorui-bi/Ariadne.git
cd Ariadne
conda env create -f environment.yml
conda activate ariadne
pip install -e .
```

最小安装：

```bash
pip install -e .
```

启用 ESM2 CeeSs 打分：

```bash
pip install -e '.[esm]'
```

## 快速开始

准备参考 FASTA 目录后运行：

```bash
ariadne run \
  --protein-folder my_proteins/ \
  --reference-dir reference_fastas/ \
  --ceess-xlsx TPS.xlsx \
  --output-dir results/
```

转录组模式会先用 Pyrodigal 预测 ORF，再进入发现阶段：

```bash
ariadne run \
  --transcriptomes sample1.fasta sample2.fasta \
  --reference-dir reference_fastas/ \
  --output-dir results_from_transcriptomes/
```

核心输出结构：

```text
results/
├── 01_discovery/
├── 02_filtering/
├── 03_classification/
└── pipeline_summary.tsv
```

湿实验验证是下游实验步骤。Ariadne 提供计算证据，并不执行酶活测定。

## 文档

- [Read the Docs](https://ariadne-platform.readthedocs.io/en/latest/)
- [Getting Started](./docs/getting-started.md)
- [Method](./docs/method.md)
- [CeeSs Classifier](./docs/esm-type.md)
- [Outputs](./docs/outputs.md)
- [CLI Reference](./docs/cli-reference.md)
- [Citation](./docs/citation.md)

## 开发

```bash
pip install -e '.[dev]'
python -m ariadne --help
mkdocs build --strict
```

## 引用

如果 Ariadne 对你的工作有帮助，请引用：

Wang, Y.; Yu, M.; Jiang, Z.; Zhou, C.; Feng, W.; Yu, K.; Ju, J.; Li, F. Identification of Cembrene Synthases Uncovers Polyphyletic Origins of 14-Membered Carbocyclic Cembranoids Biosynthesis in Corals. *J. Nat. Prod.* **2026**. DOI: [10.1021/acs.jnatprod.6c00686](https://doi.org/10.1021/acs.jnatprod.6c00686)

```bibtex
@article{wang2026cembrene,
  author    = {Wang, Yinhao and Yu, Mengmeng and Jiang, Zhaorui and Zhou, Chengyu and Feng, Wei and Yu, Ke and Ju, Jianhua and Li, Feng},
  title     = {Identification of Cembrene Synthases Uncovers Polyphyletic Origins of 14-Membered Carbocyclic Cembranoids Biosynthesis in Corals},
  journal   = {Journal of Natural Products},
  year      = {2026},
  month     = {sep},
  publisher = {American Chemical Society},
  doi       = {10.1021/acs.jnatprod.6c00686},
  url       = {https://doi.org/10.1021/acs.jnatprod.6c00686},
  note      = {Published online September 8, 2026}
}
```

## 许可证

基于 [MIT License](./LICENSE) 发布。
