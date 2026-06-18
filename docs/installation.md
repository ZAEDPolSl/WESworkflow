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
- GLnexus
- Beagle
- ANNOVAR

External resources are specified in our [**external resources guide**](docs/external_resources.md)

## 5. Run the installation check

After editing the local configuration file, run:

    bash scripts/check_installation.sh config/local_config.yaml

For the example configuration file, the script will report warnings for placeholder paths:

    bash scripts/check_installation.sh config/example_config.yaml

These warnings are expected and indicate that the user has not yet provided local paths.

The installation check verifies:

- command-line tools available in the Conda environment
- Python package imports
- R package imports
- YAML configuration syntax
- existence of configured external tools, files, and directories

This check is intended as a lightweight dry run of the installation. It does not perform full WES processing.
