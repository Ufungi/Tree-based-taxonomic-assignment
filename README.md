This pipeline performs phylogenetic placement of query sequences onto a backbone tree, refines the tree structure via divide-and-conquer methods, clusters tips into phylotypes with monophyly verification, and assigns taxonomy using a custom tree-walking algorithm to provide accurate taxonomic annotation of input sequences.

# 0. Prerequisites
* Linux environment (tested on Ubuntu 20.04)
* [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) environment

# 1. Installation
Clone **this repository** and **modified `udance`** repository.

Set conda environment
```bash
git clone https://github.com/Ufungi/Tree-based-taxonomic-assignment
cd Tree-based-taxonomic-assignment
conda env create -f environment.yaml
conda activate TM-pipeline
Rscript -e "install.packages(c('MonoPhy', 'treestats'), repos='https://cloud.r-project.org')"
```

Install modified uDance conda environment
```
git clone https://github.com/Ufungi/uDance.git
cd uDance
bash -l install.sh
```
---
### All processes below should be done in "Tree-based-taxonomic-assignment" directory
# 2. Input files


This pipeline requires three input files:

- **query.fasta**: Query sequences to be placed and classified
- **prefix_DB.fasta**: Reference backbone sequences
- **taxonomy.txt**: Taxonomic information corresponding to reference backbone sequences

Make sure paths to these files are correctly set in config.yaml.

The `Agaricales/` folder contains a minimal example input.  
You can test the pipeline using its configuration and marker files.
This example reproduces the results of Yoo et al. (2022), where genus-level classifications were fully consistent and all but one species-level classification matched the original study.
https://www.tandfonline.com/doi/full/10.1080/12298093.2022.2097364

The `Bacteria/` and `Fungi/` folders contain bacterial and fungal sequences used in the study `Cross-study discovery of functionally distinct core microbes of Tricholoma matsutake`

---

# 3. Configuration

Edit `config.yaml` to set:
- organism: Target organism or group name
- db_prefix: Prefix for input database FASTA
- threads: Number of threads to use
- chunk_size: Number of sequences per chunk
- outgroup: Comma-separated list of outgroup sequence IDs
- wdr: Working directory path
- witch, udance, fastroot: Paths to external tools

This file controls all behavior of the pipeline.

---

# 4. Run

```bash
bash run.sh
```

---
