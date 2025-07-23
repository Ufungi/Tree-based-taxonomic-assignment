This repository contains a Snakemake-based pipeline for phylogenetic analysis of fungal markers, with a specific example dataset from the order Agaricales.

📁 Repository Structure
bash
복사
편집
├── .snakemake/         # Snakemake internal files
├── Agaricales/         # Example dataset for testing
├── scripts/            # Contains custom scripts including udance modifications
├── LICENSE             # MIT license
├── config.yaml         # User-editable parameters
├── environment.yaml    # Conda environment specification
└── run.sh              # Main launcher script
⚙️ Prerequisites
conda (recommended: Miniconda)

Snakemake (installed via environment.yaml)

udance (forked and modified version – must use your fork, not the original)

🔧 Setup Instructions
1. Clone the repository
bash
복사
편집
git clone https://github.com/<your-username>/your-forked-repo.git
cd your-forked-repo
2. Create the conda environment
bash
복사
편집
conda env create -f environment.yaml
conda activate phylo_env  # or whatever name is specified in environment.yaml
3. Install modified udance (do not use the original repo)
This pipeline requires a modified version of udance. Make sure your fork includes changes to:

udance.smk

process_a_marker.sh

Clone and install your modified version:

bash
복사
편집
git clone https://github.com/<your-username>/udance.git
cd udance
# Follow installation instructions in udance repo (e.g., pip install -e . or setup.sh)
✏️ Configuration
Edit config.yaml to suit your dataset and parameters.

Key fields include:

yaml
복사
편집
markers:
  - ITS
  - LSU
threads: 8
output_dir: results/
...
Make sure paths and marker names match your actual files.

▶️ Running the Pipeline
Once everything is configured:

bash
복사
편집
bash run.sh
run.sh automatically triggers the Snakemake pipeline based on config.yaml.

🧪 Example Dataset
The Agaricales/ directory provides an example input dataset you can use to test the pipeline.

To try it:

bash
복사
편집
# (Optional) Copy config
cp Agaricales/config.yaml config.yaml

# Then run
bash run.sh
❗ Notes
run.sh relies entirely on config.yaml, so verify all parameter paths are correct.

Do not install udance from the original repository. It lacks the necessary modifications.

This pipeline is currently optimized for fungal phylogenetic markers (e.g., ITS, LSU), but can be extended to others.

📜 License
This project is licensed under the MIT License. See LICENSE for details.
