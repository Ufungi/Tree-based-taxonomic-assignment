# 1. Requirements
conda

A modified fork of udance

You must use the forked version that includes changes to udance.smk and process_a_marker.sh.

git clone https://github.com/<your-username>/udance.git
cd udance
Follow installation instructions in udance repo (e.g., pip install -e . or setup.sh)

# 2. Setup
Clone this repository and your modified udance repository.

Create and activate the conda environment:

**conda env create -f environment.yaml
conda activate phylo_env
Install your modified udance manually. Do not install the original udance.**

# 3. Configuration
Edit config.yaml to set your markers, number of threads, and input/output paths.
This file controls all behavior of the pipeline and is required before running.

# 4. Run
Use the provided shell script to launch the workflow:

**bash run.sh
This script reads from config.yaml and starts the Snakemake pipeline accordingly.
**

# 5. Example Dataset
The Agaricales/ folder contains a minimal example input.
You can test the pipeline by using its configuration and marker files.

Notes
run.sh is tightly coupled to config.yaml. Make sure to review all paths and options.

Only use your modified udance. Installing the original one will result in errors.

