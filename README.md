# 1. Requirements

- `conda`  
- `udance` (modified version)

You **must use** the forked version that includes changes to:
- `udance.smk`
- `process_a_marker.sh`

```bash
git clone https://github.com/<your-username>/udance.git
cd udance
# Follow the installation instructions in the udance repo:
# e.g., 
pip install -e .
# or
./setup.sh
```

---

# 2. Setup

Clone **this repository** and your **modified `udance`** repository.

Create and activate the conda environment:

```bash
conda env create -f environment.yaml
conda activate phylo_env
Rscript -e "install.packages(c('MonoPhy', 'treestats') repos='https://cloud.r-project.org')"
```

Then install your **modified `udance` manually**.  
Do **not** install the original `udance`.

---

# 3. Configuration

Edit `config.yaml` to set:
- your markers
- number of threads
- input/output paths

This file **controls all behavior** of the pipeline and **must be configured** before running.

---

# 4. Run

Use the provided shell script to launch the workflow:

```bash
bash run.sh
```

This script reads from `config.yaml` and starts the **Snakemake** pipeline accordingly.

---

# 5. Example Dataset

The `Agaricales/` folder contains a **minimal example input**.  
You can test the pipeline using its configuration and marker files.

---

# Notes

- `run.sh` is **tightly coupled** to `config.yaml`.  
  ✅ Be sure to review all paths and options.

- ❗ Only use your **modified** `udance`.  
  Installing the original version will **result in errors**.
