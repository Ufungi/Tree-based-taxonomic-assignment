# Tree-based Taxonomic Assignment Pipeline

## Installation

### 1. Clone this repository

``` bash
git clone https://github.com/Ufungi/Tree-based-taxonomic-assignment
cd Tree-based-taxonomic-assignment
```

### 2. Create the main conda environment

``` bash
conda env create -f environment.yaml
conda activate TM-pipeline
Rscript -e "install.packages(c('MonoPhy', 'treestats'), repos='https://cloud.r-project.org')"
```

### 3. Install the modified uDance (required for tree refinement)

``` bash
git clone https://github.com/Ufungi/uDance.git
cd uDance
bash -l install.sh
cd ..
```

------------------------------------------------------------------------

## Input Files

The pipeline requires three input files placed in `<organism>/in/`:

  File                  Description
  --------------------- ----------------------------------------------
  query.fasta           Query sequences to classify
  `<prefix>_DB.fasta`   Reference backbone sequences
  taxonomy.txt          Tab-separated: sequence ID → taxonomic ranks

**Notes:** - Sequence IDs in `<prefix>_DB.fasta` and `taxonomy.txt` must
match exactly. - Avoid special characters in filenames (underscores are
also disallowed).

### taxonomy.txt format (tab-separated)

    SEQ_ID    Kingdom    Phylum    Class    Order    Family    Genus    Species

------------------------------------------------------------------------

## Usage

### Run all steps

``` bash
conda activate TM-pipeline
bash run.sh
```

### Run a specific step

``` bash
bash run.sh 3
bash run.sh 3-6
bash run.sh --help
```

------------------------------------------------------------------------

## Example Data

The `Agaricales/` folder provides a minimal working example.

``` bash
bash run.sh
```

------------------------------------------------------------------------

## Citation

Yoo, S.N. et al. (2025). *Cross-study discovery of functionally distinct
core microbes of Tricholoma matsutake*.\
Microbiome, 13, 107.\
https://doi.org/10.1186/s40168-025-02329-x
