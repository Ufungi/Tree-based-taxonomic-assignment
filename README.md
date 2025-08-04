This pipeline performs phylogenetic placement of query sequences onto a backbone tree, refines the tree structure via divide-and-conquer methods, clusters tips into phylotypes with monophyly verification, and assigns taxonomy using a custom tree-walking algorithm to provide accurate taxonomic annotation of input sequences.

# 1. Installation
Clone **this repository** and **modified `udance`** repository.

Set conda environment
```bash
git clone https://github.com/Ufungi/Tree-based-taxonomic-assignment
conda env create -f environment.yaml
conda activate phylo_env
Rscript -e "install.packages(c('MonoPhy', 'treestats'), repos='https://cloud.r-project.org')"
```

Install modified uDance conda environment
```
git clone https://github.com/Ufungi/uDance.git
cd uDance
bash install.sh
```
---

# 2. Input files

This pipeline requires three input files:

- **query.fasta**: Query sequences to be placed and classified
- **prefix_DB.fast**a: Reference backbone sequences for phylogenetic placement
- **taxonomy.txt**: Taxonomic information corresponding to reference sequences

Make sure paths to these files are correctly set in config.yaml.

The `Agaricales/` folder contains a minimal example input.  
You can test the pipeline using its configuration and marker files.

---

# 3. Configuration

Edit `config.yaml` to set:
- marker genes
- number of threads
- input/output paths

This file controls all behavior of the pipeline.

---

# 4. Run

```bash
bash run.sh
```

This script reads from `config.yaml` and starts the pipeline.

---
