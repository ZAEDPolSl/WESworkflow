#!/usr/bin/env bash
set -euo pipefail

# This script computes gene-level features per chromosome:
# (1) takes functionally filtered and CADD-annotated original VCF files,
# (2) combines them with corresponding filtered imputed VCF files,
# (3) computes per-sample, per-gene features using cal_features_multi.py.
# Usage: bash Variant_to_gene/gene_aggregation.sh config/local_config.yaml

CONFIG="${1:-config/local_config.yaml}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"

if [[ "$CONFIG" != /* ]]; then
	CONFIG="$REPO_DIR/$CONFIG"
fi

READ_CONFIG="$REPO_DIR/scripts/read_config.py"

resolve_path() {
	local path="$1"

	if [[ "$path" == /* ]]; then
		echo "$path"
	else
		echo "$REPO_DIR/$path"
	fi
}

if [[ ! -f "$CONFIG" ]]; then
	echo "ERROR: Config file not found: $CONFIG"
	exit 1
fi

RESULTS_DIR="$(python "$READ_CONFIG" "$CONFIG" directories.results_dir)"
PARALLEL_JOBS="$(python "$READ_CONFIG" "$CONFIG" parameters.feature_parallel_jobs)"

RESULTS_DIR="$(resolve_path "$RESULTS_DIR")"

IMPUTED_DIR="$RESULTS_DIR/Imputation/Output_Chromosomes/Beagle_1kG_imputed"
ORIGINAL_DIR="$RESULTS_DIR/Annotation/CADD_scored"
FEATURES_DIR="$RESULTS_DIR/Features"
FEATURE_SCRIPT="$REPO_DIR/Variant_to_gene/cal_features_multi.py"

CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

if [[ ! -d "$IMPUTED_DIR" ]]; then
	echo "ERROR: Imputed VCF directory not found: $IMPUTED_DIR"
	exit 1
fi

if [[ ! -d "$ORIGINAL_DIR" ]]; then
	echo "ERROR: CADD-annotated VCF directory not found: $ORIGINAL_DIR"
	exit 1
fi

if [[ ! -f "$FEATURE_SCRIPT" ]]; then
	echo "ERROR: Feature calculation script not found: $FEATURE_SCRIPT"
	exit 1
fi

mkdir -p "$FEATURES_DIR"

calculate_features() {
	local chrom="$1"
	local original_vcf="$ORIGINAL_DIR/chr${chrom}.hg38_multianno.filtered_CADD_scored.vcf.gz"
	local imputed_vcf="$IMPUTED_DIR/chr${chrom}.filtered.vcf.gz"
	local out_dir="$FEATURES_DIR/chr${chrom}"

	echo "Calculating features for chr${chrom}"

	if [[ ! -f "$original_vcf" ]]; then
		echo "ERROR: Original CADD-annotated VCF not found: $original_vcf"
		exit 1
	fi

	if [[ ! -f "$imputed_vcf" ]]; then
		echo "ERROR: Filtered imputed VCF not found: $imputed_vcf"
		exit 1
	fi

	mkdir -p "$out_dir"

	python "$FEATURE_SCRIPT" \
		--original "$original_vcf" \
		--imputed "$imputed_vcf" \
		--out_dir "$out_dir"
}

export -f calculate_features
export IMPUTED_DIR ORIGINAL_DIR FEATURES_DIR FEATURE_SCRIPT

printf "%s\n" "${CHROMOSOMES[@]}" | \
	parallel --halt soon,fail=1 -j "$PARALLEL_JOBS" \
	bash -c 'calculate_features "$1"' _ {}

echo "Gene-level feature calculation completed."
