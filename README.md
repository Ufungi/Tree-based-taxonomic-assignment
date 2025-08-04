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

# 2. Configuration

Edit `config.yaml` to set:
- marker genes
- number of threads
- input/output paths

This file controls all behavior of the pipeline.

---

# 3. Run

```bash
bash run.sh
```

This script reads from `config.yaml` and starts the pipeline.

---

# 4. Example Dataset

The `Agaricales/` folder contains a minimal example input.  
You can test the pipeline using its configuration and marker files.
