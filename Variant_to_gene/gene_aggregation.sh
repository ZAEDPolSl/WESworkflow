#!/bin/bash

# This script computes gene-level features per chromosome:
# (1) Takes functionally filtered and CADD-annotated VCF (original)
# (2) Combines it with the corresponding imputed VCF
# (3) Computes per-sample, per-gene features using cal_features_multi.py
#
# Input: per-chromosome original (CADD-annotated) and imputed VCF files
# Output: per-chromosome feature matrices stored in separate directories

WORKDIR="/path/to/results/dir" # all the downstream subdirectories with intermediate files will be created here
IMPUTED_DIR="$WORKDIR/Imputation/Output_Chromosomes/Beagle_1kG_imputed"
ORIGINAL_DIR="$WORKDIR/Annotation/CADD_scored"

mkdir -p $WORKDIR/Features


# step 6: calculate features ---------------------------------------------------------

export IMPUTED_DIR ORIGINAL_DIR WORKDIR
parallel -j 23 '
  mkdir -p "$WORKDIR/Features/chr{1}" &&
  cal_features_multi.py \
    --original $ORIGINAL_DIR/chr{1}.hg38_multianno.filtered_CADD_scored.vcf.gz \
    --imputed  $IMPUTED_DIR/chr{1}.filtered.vcf.gz \
    --out_dir  $WORKDIR/Features/chr{1} \
' ::: {1..22} X 