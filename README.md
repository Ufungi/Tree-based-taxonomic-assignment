# Tree-based Taxonomic Assignment Pipeline

<p align="center">
  <img src="https://img.shields.io/badge/platform-Linux-blue?logo=linux&logoColor=white" alt="Platform">
  <img src="https://img.shields.io/badge/conda-supported-green?logo=anaconda&logoColor=white" alt="Conda">
  <img src="https://img.shields.io/badge/snakemake-compatible-blue?logo=snakemake&logoColor=white" alt="Snakemake">
  <img src="https://img.shields.io/badge/language-Bash%20%7C%20R%20%7C%20Python-informational?logo=gnu-bash&logoColor=white" alt="Language">
  <img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="License">
  <img src="https://img.shields.io/badge/DOI-10.1186%2Fs40168--025--02329--x-orange" alt="DOI">
</p>

> **Accurate taxonomic annotation of metagenomic/amplicon sequences via phylogenetic placement, tree refinement, and LCA-based assignment.**

This pipeline places query sequences onto a reference backbone tree, refines the topology using divide-and-conquer methods (uDance), clusters tips into phylotypes with monophyly verification, and assigns taxonomy through a custom tree-walking LCA algorithm — providing robust, phylogeny-aware classification for bacteria, fungi, and other microbial groups.

---

## Table of Contents

- [Key Features](#-key-features)
- [Pipeline Overview](#-pipeline-overview)
- [Prerequisites](#-prerequisites)
- [Installation](#-installation)
- [Input Files](#-input-files)
- [Configuration](#-configuration)
- [Usage](#-usage)
- [Example Data](#-example-data)
- [Tools Used](#-tools-used)
- [Citation](#-citation)

---

## ✨ Key Features

- **Phylogenetic placement** — Places query sequences onto a reference tree using [App-SpaM](https://github.com/matthiasblanke/App-SpaM) (alignment-free, fast)
- **Tree refinement** — Expands the backbone tree with [uDance](https://github.com/Ufungi/uDance) divide-and-conquer subtree updates
- **Phylotype clustering** — Groups sequences by phylogenetic distance using [TreeCluster](https://github.com/niemasd/TreeCluster)
- **Long-branch detection** — Removes spurious sequences via [TreeShrink](https://github.com/uym2/TreeShrink)
- **LCA taxonomy** — Assigns taxonomy using a monophyly-verified, Last Common Ancestor algorithm in R
- **Modular execution** — Run all 7 steps or jump to any specific step or range

---

## Pipeline Overview


Input: query.fasta + reference DB (FASTA + taxonomy)
│
[Step 1] MSA (WITCH) → Tree inference (VeryFastTree) → Branch rescaling (RAxML-NG)
│
[Step 2] Phylogenetic placement of query sequences (App-SpaM / gappa chunkify)
│
[Step 3] jplace validation and cleanup
│
[Step 4] Tree refinement with uDance (divide-and-conquer subtree expansion)
│
[Step 5] Phylotype clustering (TreeCluster, threshold = 0.02)
│
[Step 6] Long-branch detection and removal (TreeShrink)
│
[Step 7] LCA-based taxonomy assignment (R: MonoPhy + treestats)
│
Output: taxonomy table with assignments per query sequence


---

## Prerequisites

- **OS**: Linux (tested on Ubuntu 20.04)
- **Conda**: [Miniconda or Anaconda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html)
- **Storage**: ~10 GB free space for tools + reference data
- **RAM**: 16 GB minimum recommended (32+ GB for large datasets)

---

## Installation

### 1. Clone this repository

```bash
git clone https://github.com/Ufungi/Tree-based-taxonomic-assignment
cd Tree-based-taxonomic-assignment

2. Create the main conda environment
conda env create -f environment.yaml
conda activate TM-pipeline
Rscript -e "install.packages(c('MonoPhy', 'treestats'), repos='https://cloud.r-project.org')"

3. Install the modified uDance (required for tree refinement)
git clone https://github.com/Ufungi/uDance.git
cd uDance
bash -l install.sh
cd ..

Input Files
The pipeline requires three input files placed in <organism>/in/:

File	Description
query.fasta	Query sequences to classify (no special characters in filename)
<prefix>_DB.fasta	Reference backbone sequences
taxonomy.txt	Tab-separated: sequence ID → taxonomic ranks
Note: Sequence IDs in <prefix>_DB.fasta and taxonomy.txt must match exactly.
Query sequence filenames must not contain special characters (underscores are also disallowed).

taxonomy.txt format (tab-separated):

SEQ_ID    Kingdom    Phylum    Class    Order    Family    Genus    Species

Configuration
Edit config.yaml before running:

# Organism/group name — used to locate input files and name outputs
organism: Agaricales

# Prefix of the reference FASTA file (without .fasta extension)
db_prefix: Agaricales_DB

# Number of CPU threads
threads: 32

# Sequences per chunk for phylogenetic placement (or 'auto')
chunk_size: 500

# Outgroup sequence IDs for RAxML-NG rooting (comma-separated)
outgroup: NR_103631, NR_178158, NR_152951

# Absolute path to this repository (cloned folder)
wdr: /absolute/path/to/Tree-based-taxonomic-assignment

# Paths to external tools
witch: /path/to/witch.py          # which witch.py
udance: /path/to/uDance           # directory of cloned uDance
fastroot: /path/to/FastRoot.py    # which FastRoot.py

Tip: Run which witch.py and which FastRoot.py after activating the conda environment to find the tool paths.

Also update the CONFIG_GLOBAL path at the top of run.sh to point to your config.yaml:

# In run.sh, line 8 — change to your absolute path:
export CONFIG_GLOBAL="/absolute/path/to/Tree-based-taxonomic-assignment/config.yaml"

Usage
Run all steps
conda activate TM-pipeline
bash run.sh

Run a specific step
bash run.sh 3        # Run only step 3 (fix_jplace_file)
bash run.sh 3-6      # Run steps 3 through 6
bash run.sh --help   # Show step descriptions

Step reference
Step	Function	Description
1	align_and_build_tree	Align reference DB with WITCH → infer tree with VeryFastTree → rescale branches with RAxML-NG
2	phylogenetic_placement	Chunkify query sequences and place onto backbone tree with App-SpaM
3	fix_jplace_file	Validate and repair jplace output (remove NaN entries, normalize format)
4	tree_refinement	Expand backbone with placed queries using uDance subtree updates
5	cluster_analysis	Root subtrees with FastRoot, cluster tips with TreeCluster (threshold = 0.02)
6	detect_long_branches	Identify and remove long-branched tips using TreeShrink
7	run_taxonomy_script	Assign taxonomy per query via LCA R script using MonoPhy monophyly checking
Each step is resumable — if output files from a step already exist, that step is skipped automatically.

Example Data
The Agaricales/ folder provides a minimal working example. It reproduces the genus- and species-level classifications from Yoo et al. (2022) (genus-level: 100% consistent; species-level: all but one matched):

# Already configured for the example — just run:
bash run.sh

Reference: Yoo et al. (2022) Mycobiology — https://www.tandfonline.com/doi/full/10.1080/12298093.2022.2097364

The Bacteria/ and Fungi/ folders contain the bacterial and fungal reference sequences used in our study of Tricholoma matsutake core microbiome (see Citation).

Tools Used
Tool	Role	Reference
WITCH	Reference multiple sequence alignment	Shen et al.
VeryFastTree	Fast maximum-likelihood tree inference	
RAxML-NG	Branch length re-estimation & rooting	Kozlov et al.
App-SpaM	Alignment-free phylogenetic placement	Blanke & Morgenstern
gappa	jplace chunkify/unchunkify utility	Czech et al.
uDance	Divide-and-conquer tree refinement	Balaban et al.
FastRoot	Minimum-variance rooting	Mai et al.
TreeCluster	Phylotype clustering	Balaban et al.
TreeShrink	Long-branch detection	Mai & Mirarab
MonoPhy (R)	Monophyly assessment	Schwery & O'Meara
Citation
If you use this pipeline in your research, please cite:

Yoo, S.N. et al. (2025). Cross-study discovery of functionally distinct core microbes of Tricholoma matsutake.
Microbiome, 13, 107.
https://doi.org/10.1186/s40168-025-02329-x

@article{yoo2025tricholoma,
  title     = {Cross-study discovery of functionally distinct core microbes of \textit{Tricholoma matsutake}},
  author    = {Yoo, S.N. and others},
  journal   = {Microbiome},
  volume    = {13},
  pages     = {107},
  year      = {2025},
  doi       = {10.1186/s40168-025-02329-x},
  url       = {https://doi.org/10.1186/s40168-025-02329-x}
}

Note: Please verify author names and volume/page numbers against the published article at https://link.springer.com/article/10.1186/s40168-025-02329-x and update the BibTeX accordingly.

<p align="center"> <sub>Keywords: taxonomic assignment · phylogenetic placement · microbiome · ITS · 16S rRNA · amplicon · metagenomics · fungi · bacteria · <em>Tricholoma matsutake</em> · LCA · tree-based classification · uDance · App-SpaM · WITCH</sub> </p> ```
