# Installation and environment setup

This document describes how to prepare the computational environment required to run WESworkflow. The workflow uses a single Conda environment for command-line tools, Python packages, and most R packages. A small number of R packages are installed with an additional R script.

## 1. Clone the repository

    git clone https://github.com/ZAEDPolSl/WESworkflow.git
    cd WESworkflow

## 2. Create the Conda environment

    conda env create -f envs/wes-env.yml
    conda activate wes-env

If the environment already exists, update it with:

    conda env update -n wes-env -f envs/wes-env.yml --prune
    conda activate wes-env

## 3. Install additional R packages

Some R packages are installed separately to avoid Conda dependency conflicts.

    R_LIBS= R_LIBS_USER= R_LIBS_SITE= R_PROFILE_USER=/dev/null Rscript --vanilla scripts/install_R_packages.R

The explicit R_LIBS settings are used to prevent R from loading packages from a user-level library outside the Conda environment.

## 4. Configure external tools and resources

Copy the example configuration file:

    cp config/example_config.yaml config/local_config.yaml

Edit config/local_config.yaml and replace all /path/to/... placeholders with local paths to the required tools and resources.

The following external tools are expected to be configured manually:

- DeepVariant
- Beagle
- ANNOVAR

External resources are specified in our [**external resources guide**](external_resources.md)

## 5. Computational resource configuration

The default runtime parameters in `config/example_config.yaml` were configured for a system with **48 physical CPU cores (96 threads), 503 GiB RAM, and two NVIDIA RTX 6000 Ada GPUs (49 GB VRAM each)**. Users should adjust these settings according to their available hardware.

Parallelization is controlled by two types of parameters:
- `threads` – CPU threads used within an individual process.
- `jobs` – independent processes executed concurrently across samples or chromosomes.

When configuring parallel execution, consider the combined CPU and memory requirements (`jobs × threads` and `jobs × memory per process`), as well as storage throughput.

Computational requirements vary across workflow stages. BWA benefits from multithreading, while FastQC parallelization may be limited by disk throughput. DeepVariant uses CPU threads for read processing and optionally GPU acceleration for variant calling. Beagle is both CPU- and memory-intensive, with genotype imputation parallelized across chromosomes. Annotation and gene-level aggregation also involve substantial disk I/O.

Chromosome-level parallelization is limited to 23 jobs (chromosomes 1–22 and X), so increasing the corresponding `jobs` parameters beyond this number provides no additional benefit.

## 6. Run the installation check

After editing the local configuration file, run:

```bash
bash scripts/check_installation.sh config/local_config.yaml
```

For the example configuration file, the script will report warnings for placeholder paths:

```bash
bash scripts/check_installation.sh config/example_config.yaml
```

These warnings are expected and indicate that the user has not yet provided local paths.

The installation check verifies:

- command-line tools available in the Conda environment
- Python package imports
- R package imports
- YAML configuration syntax
- existence of configured external tools, files, and directories

This check is intended as a lightweight dry run of the installation.
