# 🧬 Ariadne: Unraveling the Labyrinth of Animal Terpene Synthases

<div align="center">
  <img src="./fig/logo.png" alt="Ariadne Pipeline Logo" width="300">
</div>

## 📖 About The Project

Ariadne is a comprehensive bioinformatics pipeline designed to mine, filter, phylogenetically classify, and functionally annotate animal-derived Terpene Synthases (TPS). Our primary focus is the identification of novel cembrene-producing enzymes from corals, insects, and other eukaryotic non-plant lineages.

In Greek mythology, Ariadne provided the thread that guided Theseus safely out of the Minotaur's labyrinth. Similarly, the Ariadne pipeline utilizes targeted Hidden Markov Models (HMM) as its "golden thread," navigating the immense complexity of raw transcriptome data (represented by the outer ring of sequence data in the logo) to pinpoint highly conserved functional motifs.

### 🎨 About The Logo Design

The Ariadne logo is a circular badge designed to visually represent the entire bioinformatics pipeline and its core mission:

- **The Labyrinth (Outer Ring):** A deep blue ring with intricate DNA and protein sequence codes (e.g., `GEEICK...`) representing the complex and often obscure raw transcriptome data.
- **The Golden Thread (Central Path):** A luminous golden line that originates from a DNA double helix and traces a path through a series of peptide residues (`A S C V G I G G G S T G P ...`), representing the HMM-driven search through sequence space.
- **The Minotaur's Lair (Core):** The thread culminates in the core of the labyrinth, where the central target—the highly conserved and critical catalytic motif **CFDVLAF**—is revealed in a distinct, focused circle. This is the precise location of functional annotation.
- **Biological Context:** Faint silhouettes of a coral fan and insect wings frame the labyrinth, underscoring our focus on animal-derived TPS.

## ⚙️ Pipeline Architecture

The Ariadne pipeline is divided into four main phases:

### Phase 1: Labyrinth Entry (Translation & HMM Discovery)

Transforms raw transcriptome assemblies into translated protein space and performs homology-based sequence mining.

- **ORF Prediction:** Utilizes Prodigal to translate the whole transcriptome data into amino acid sequences (output routed to the protein/ directory).
- **HMM Searching:** Scans the proteome using customized TPS Hidden Markov Models to generate candidate hit files (xx.out).
- **Sequence Extraction:** Automates the extraction of full candidate sequences using the custom perl script (`perl /data/yym/fasta_extract_mult.pl`).

### Phase 2: Refining the Thread (Quality Control & De-replication)

Filters the raw candidates to retain only high-quality, non-redundant, and structurally complete sequences.

- **De-duplication:** Runs `DupRemover.py` to collapse redundant sequences, retaining only unique isoforms (≥ 95% similarity threshold).
- **Coverage Filtering:** Uses `filterxx.py` to strictly remove truncated sequences or those with low HMM coverage (cov < 10).
- _Note: Manual curation is highly recommended at this stage to ensure the retained sequences possess complete N- and C-termini._

### Phase 3: Mapping the Maze (Phylogenetic & Structural Profiling)

Determines the evolutionary origin of the refined candidates to confidently classify them as animal-derived TPS.

- **Phylogenetic Contextualization:** Constructs phylogenetic trees aligning candidates against known TPS databases including bacteria, fungi, plants, insects (insect-ggpps-fpps), and corals (coralTPS-cembrene).
- **3D Landscape Visualization:** Projects sequences into a 3D functional space (inspired by AFLP-finfer methodology) to visually confirm clustering within specific taxonomic clades.

### Phase 4: The Minotaur's Lair (Motif-Based Functional Annotation)

Pinpoints specific catalytic capabilities, identifying true cembrene-producing synthases.

- **Micro-Clustering:** Performs high-resolution sub-clustering based on the highly conserved ~210 aa region, specifically targeting the **CFDVLAF** motif (derived from S_CdTC-2 and prominently featured in the Logo's core).
- **Motif Visualization:** Generates motif block diagrams mapping known cembrene producers against novel sequences.
- _Note: This step utilizes the MEME Suite. Due to dependency requirements, MEME is executed in a dedicated, isolated Conda environment._

## 🛠️ Prerequisites & Installation

To run Ariadne successfully, ensure you have the following software and dependencies installed:

- **Core Tools:** Prodigal, HMMER3
- **Python Dependencies:** Biopython, pandas, numpy (for custom Python scripts)
- **Perl:** Base installation (for sequence extraction scripts)
- **MEME Suite:** Requires a separate Conda environment due to legacy dependencies.

## 🚀 Quick Start / Usage

(Here you can add specific command-line examples for users. Placeholder below)

## 📄 License

Distributed under the MIT License. See `LICENSE` for more information.
