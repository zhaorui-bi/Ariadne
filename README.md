# 🧬 Ariadne: Unraveling the Labyrinth of Animal Terpene Synthases

## 📖 About The Project
In Greek mythology, Ariadne provided the thread that guided Theseus safely out of the Minotaur's labyrinth. Similarly, the Ariadne Pipeline utilizes targeted Hidden Markov Models (HMM) and highly conserved catalytic motifs (such as CFDVLAF) as the "golden thread" to navigate the immense complexity of raw transcriptome data.

This comprehensive bioinformatics workflow is specifically designed to mine, filter, phylogenetically classify, and functionally annotate animal-derived Terpene Synthases (TPS)—with a special focus on identifying novel cembrene-producing enzymes from corals, insects, and other eukaryotic non-plant lineages.

## ⚙️ Pipeline Architecture
The Ariadne pipeline is divided into four main phases:

### Phase 1: Labyrinth Entry (Translation & HMM Discovery)
Transforms raw transcriptome assemblies into translated protein space and performs homology-based sequence mining.

ORF Prediction: Utilizes Prodigal to translate the whole transcriptome data into amino acid sequences (output routed to the protein/ directory).

HMM Searching: Scans the proteome using customized TPS Hidden Markov Models to generate candidate hit files (xx.out).

Sequence Extraction: Automates the extraction of full candidate sequences using the custom perl script (perl /data/yym/fasta_extract_mult.pl).

### Phase 2: Refining the Thread (Quality Control & De-replication)
Filters the raw candidates to retain only high-quality, non-redundant, and structurally complete sequences.

De-duplication: Runs DupRemover.py to collapse redundant sequences, retaining only unique isoforms (≥ 95% similarity threshold).

Coverage Filtering: Uses filterxx.py to strictly remove truncated sequences or those with low HMM coverage (cov < 10).

Note: Manual curation is highly recommended at this stage to ensure the retained sequences possess complete N- and C-termini.

### Phase 3: Mapping the Maze (Phylogenetic & Structural Profiling)
Determines the evolutionary origin of the refined candidates to confidently classify them as animal-derived TPS.

Phylogenetic Contextualization: Constructs phylogenetic trees aligning candidates against known TPS databases including bacteria, fungi, plants, insects (insect-ggpps-fpps), and corals (coralTPS-cembrene).

3D Landscape Visualization: Projects sequences into a 3D functional space (inspired by AFLP-finfer methodology) to visually confirm clustering within specific taxonomic clades.

### Phase 4: The Minotaur's Lair (Motif-Based Functional Annotation)
Pinpoints specific catalytic capabilities, identifying true cembrene-producing synthases.

Micro-Clustering: Performs high-resolution sub-clustering based on the highly conserved ~210 aa region, specifically targeting the CFDVLAF motif (derived from S_CdTC-2).

Motif Visualization: Generates motif block diagrams mapping known cembrene producers against novel sequences.

Note: This step utilizes the MEME Suite. Due to dependency requirements, MEME is executed in a dedicated, isolated Conda environment.

## 🛠️ Prerequisites & Installation
To run Ariadne successfully, ensure you have the following software and dependencies installed:

Core Tools: Prodigal, HMMER3

Python Dependencies: Biopython, pandas, numpy (for custom Python scripts)

Perl: Base installation (for sequence extraction scripts)

MEME Suite: Requires a separate Conda environment due to legacy dependencies.

## 🚀 Quick Start / Usage
(Here you can add specific command-line examples for users. Placeholder below)

## 📄 License
Distributed under the MIT License. See LICENSE for more information.
