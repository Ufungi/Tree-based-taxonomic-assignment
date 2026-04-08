# Phylogenetic Tree-based Taxonomic Assignment of Heterogenous Amplicon Sequences

<p align="center">
  <img src="https://img.shields.io/badge/platform-Linux-blue?logo=linux&logoColor=white" alt="Platform">
  <img src="https://img.shields.io/badge/conda-supported-green?logo=anaconda&logoColor=white" alt="Conda">
  <img src="https://img.shields.io/badge/snakemake-compatible-blue?logo=snakemake&logoColor=white" alt="Snakemake">
  <img src="https://img.shields.io/badge/language-Bash%20%7C%20R%20%7C%20Python-informational?logo=gnu-bash&logoColor=white" alt="Language">
  <img src="https://img.shields.io/badge/license-MIT-lightgrey" alt="License">
  <img src="https://img.shields.io/badge/DOI-10.1186%2Fs40168--025--02329--x-orange" alt="DOI">
</p>

> **Accurate taxonomic annotation of amplicon sequences via phylogenetic placement, tree refinement, and LCA-based assignment.**

This pipeline places query sequences onto a reference backbone tree, refines the topology using divide-and-conquer methods (uDance), clusters tips into phylotypes with monophyly verification, and assigns taxonomy through a custom tree-walking LCA algorithm — providing robust, phylogeny-aware classification for bacteria, fungi, and other microbial groups.

## ✨ Key Features

- **Phylogenetic placement** — Places query sequences onto a reference tree using [App-SpaM](https://github.com/matthiasblanke/App-SpaM) (alignment-free, fast)
- **Tree refinement** — Expands the backbone tree with [uDance](https://github.com/Ufungi/uDance) divide-and-conquer subtree updates
- **Phylotype clustering** — Groups sequences by phylogenetic distance using [TreeCluster](https://github.com/niemasd/TreeCluster)
- **Long-branch detection** — Removes spurious sequences via [TreeShrink](https://github.com/uym2/TreeShrink)
- **LCA taxonomy** — Assigns taxonomy using a monophyly-verified, Last Common Ancestor algorithm in R
- **Modular execution** — Run all 7 steps or jump to any specific step or range

---

## Pipeline Overview
<img width="769" height="703" alt="Workflow_github" src="https://github.com/user-attachments/assets/5c6ca91d-2283-4907-87f6-890fbea65d04" />

```
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
```

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
```

### 2. Create the main conda environment

```bash
conda env create -f environment.yaml
conda activate TM-pipeline
Rscript -e "install.packages(c('MonoPhy', 'treestats'), repos='https://cloud.r-project.org')"
```

### 3. Install the modified uDance (required for tree refinement)

```bash
git clone https://github.com/Ufungi/uDance.git
cd uDance
bash -l install.sh
cd ..
```

---

## Input Files

The pipeline requires three input files placed in `<organism>/in/`:

| File | Description |
|------|-------------|
| `query.fasta` | Query sequences to classify (no special characters in filename) |
| `<prefix>_DB.fasta` | Reference backbone sequences |
| `taxonomy.txt` | Tab-separated: sequence ID → taxonomic ranks |

**Note:** Sequence IDs in `<prefix>_DB.fasta` and `taxonomy.txt` must match exactly.  
Query sequence filenames must not contain special characters (underscores are also disallowed).

`taxonomy.txt` format (tab-separated):

```
SEQ_ID    Kingdom    Phylum    Class    Order    Family    Genus    Species
```

---

## Configuration

Edit `config.yaml` before running:

```yaml
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
```

**Tip:** Run `which witch.py` and `which FastRoot.py` after activating the conda environment to find the tool paths.

Also update the `CONFIG_GLOBAL` path at the top of `run.sh` to point to your `config.yaml`:

```bash
# In run.sh, line 8 — change to your absolute path:
export CONFIG_GLOBAL="/absolute/path/to/Tree-based-taxonomic-assignment/config.yaml"
```

---

## Usage

### Run all steps

```bash
conda activate TM-pipeline
bash run.sh
```

### Run a specific step

```bash
bash run.sh 3        # Run only step 3 (fix_jplace_file)
bash run.sh 3-6      # Run steps 3 through 6
bash run.sh --help   # Show step descriptions
```

### Step reference

| Step | Function | Description |
|------|----------|-------------|
| 1 | `align_and_build_tree` | Align reference DB with WITCH → infer tree with VeryFastTree → rescale branches with RAxML-NG |
| 2 | `phylogenetic_placement` | Chunkify query sequences and place onto backbone tree with App-SpaM |
| 3 | `fix_jplace_file` | Validate and repair jplace output (remove NaN entries, normalize format) |
| 4 | `tree_refinement` | Expand backbone with placed queries using uDance subtree updates |
| 5 | `cluster_analysis` | Root subtrees with FastRoot, cluster tips with TreeCluster (threshold = 0.02) |
| 6 | `detect_long_branches` | Identify and remove long-branched tips using TreeShrink |
| 7 | `run_taxonomy_script` | Assign taxonomy per query via LCA R script using MonoPhy monophyly checking |

Each step is resumable — if output files from a step already exist, that step is skipped automatically.

---

## Example Data

The `Agaricales/` folder provides a minimal working example. It reproduces the genus- and species-level classifications from Yoo et al. (2022) (genus-level: 100% consistent; species-level: all but one matched):

```bash
# Already configured for the example — just run:
bash run.sh
```

Reference: Yoo et al. (2022) *Mycobiology* — https://www.tandfonline.com/doi/full/10.1080/12298093.2022.2097364

The `Bacteria/` and `Fungi/` folders contain the bacterial and fungal reference sequences used in our study of *Tricholoma matsutake* core microbiome (see Citation).

---

## Tools Used

| Tool | Role | Reference |
|------|------|-----------|
| WITCH | Reference multiple sequence alignment | Shen et al. (2022) *J. Comput. Biol.* 29:782–801 |
| VeryFastTree | Fast maximum-likelihood tree inference | Piñeiro et al. (2020) *Bioinformatics* 36:4658–4659 |
| RAxML-NG | Branch length re-estimation & rooting | Kozlov et al. (2019) *Bioinformatics* 35:4453–4455 |
| App-SpaM | Alignment-free phylogenetic placement | Blanke & Morgenstern (2021) *Bioinform. Adv.* 1:vbab027 |
| gappa | jplace chunkify/unchunkify utility | Czech et al. (2020) *Bioinformatics* 36:3263–3265 |
| uDance | Divide-and-conquer tree refinement | Balaban et al. (2024) *Nat. Biotechnol.* 42:768–777 |
| FastRoot | Minimum-variance rooting | Mai et al. (2017) *PLOS ONE* 12:e0182238 |
| TreeCluster | Phylotype clustering | Balaban et al. (2019) *PLOS ONE* 14:e0221068 |
| TreeShrink | Long-branch detection | Mai & Mirarab (2018) *BMC Genomics* 19:272 |
| MonoPhy (R) | Monophyly assessment | Schwery & O'Meara (2016) *PeerJ Comput. Sci.* 2:e56 |

---

## Citation

If you use this pipeline in your research, please cite:

**Pipeline:**

Yoo, S., Seo, C.W. & Lim, Y.W. (2026). Functionally distinct core microbes of *Tricholoma matsutake* revealed by cross-study analysis. *Microbiome*, 14, 58. https://doi.org/10.1186/s40168-025-02329-x

**Tools:**

Shen, C., Park, M. & Warnow, T. (2022). WITCH: improved multiple sequence alignment through weighted consensus hidden Markov model alignment. *J. Comput. Biol.*, 29, 782–801. https://doi.org/10.1089/cmb.2021.0585

Piñeiro, C., Abuín, J.M. & Pichel, J.C. (2020). VeryFastTree: speeding up the estimation of phylogenies for large alignments through parallelization and vectorization strategies. *Bioinformatics*, 36, 4658–4659. https://doi.org/10.1093/bioinformatics/btaa582

Kozlov, A.M., Darriba, D., Flouri, T., Morel, B. & Stamatakis, A. (2019). RAxML-NG: a fast, scalable and user-friendly tool for maximum likelihood phylogenetic inference. *Bioinformatics*, 35, 4453–4455. https://doi.org/10.1093/bioinformatics/btz305

Blanke, M. & Morgenstern, B. (2021). App-SpaM: phylogenetic placement of short reads without sequence alignment. *Bioinform. Adv.*, 1, vbab027. https://doi.org/10.1093/bioadv/vbab027

Czech, L., Barbera, P. & Stamatakis, A. (2020). Genesis and Gappa: processing, analyzing and visualizing phylogenetic (placement) data. *Bioinformatics*, 36, 3263–3265. https://doi.org/10.1093/bioinformatics/btaa070

Balaban, M., Jiang, Y., Zhu, Q., McDonald, D., Knight, R. & Mirarab, S. (2024). Generation of accurate, expandable phylogenomic trees with uDance. *Nat. Biotechnol.*, 42, 768–777. https://doi.org/10.1038/s41587-023-01868-8

Mai, U., Sayyari, E. & Mirarab, S. (2017). Minimum variance rooting of phylogenetic trees and implications for species tree reconstruction. *PLOS ONE*, 12, e0182238. https://doi.org/10.1371/journal.pone.0182238

Balaban, M., Moshiri, N., Mai, U., Jia, X. & Mirarab, S. (2019). TreeCluster: clustering biological sequences using phylogenetic trees. *PLOS ONE*, 14, e0221068. https://doi.org/10.1371/journal.pone.0221068

Mai, U. & Mirarab, S. (2018). TreeShrink: fast and accurate detection of outlier long branches in collections of phylogenetic trees. *BMC Genomics*, 19, 272. https://doi.org/10.1186/s12864-018-4620-2

Schwery, O. & O'Meara, B.C. (2016). MonoPhy: a simple R package to find and visualize monophyly issues. *PeerJ Comput. Sci.*, 2, e56. https://doi.org/10.7717/peerj-cs.56

---

<p align="center">
  <sub>Keywords: taxonomic assignment · phylogenetic placement · microbiome · ITS · 16S rRNA · amplicon · metagenomics · fungi · bacteria · <em>Tricholoma matsutake</em> · LCA · tree-based classification · uDance · App-SpaM · WITCH</sub>
</p>
