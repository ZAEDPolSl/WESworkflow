# Example downstream data

This directory contains a minimal mock dataset for checking downstream file formats and repository setup.

The example data are not intended to reproduce full WES processing. They are small synthetic feature files used to verify that the repository structure, metadata format, and downstream feature file format are valid.

## Files

- sample_path_map_example.tsv - example metadata file with required columns:
  - Sample
  - Dataset
  - Kit

- features/ - example per-sample gene-level feature files arranged by chromosome.

Each .feature.txt file must contain at least two columns:

Gene    CADD_weighted_avg_AF

## Run the example dry run

From the repository root:

    bash scripts/run_example_downstream.sh

The script copies the example feature files into example_results/Features/ and validates:

- feature file presence
- required feature columns
- required metadata columns
- consistency between feature sample names and metadata sample names

Generated files are written to example_results/, which is ignored by Git.
