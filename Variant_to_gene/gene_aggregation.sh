#!/usr/bin/env bash
set -euo pipefail

# This script calculates gene-level features from CADD-scored original VCFs.
# If a filtered imputed VCF is available for a chromosome, it is used as an additional input.
# If the filtered imputed VCF is missing, features are calculated from the original VCF only.
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

read_cfg() {
	python "$READ_CONFIG" "$CONFIG" "$1"
}

if [[ ! -f "$CONFIG" ]]; then
	echo "ERROR: Config file not found: $CONFIG" >&2
	exit 1
fi

if [[ ! -f "$READ_CONFIG" ]]; then
	echo "ERROR: Config reader not found: $READ_CONFIG" >&2
	exit 1
fi

RESULTS_DIR="$(resolve_path "$(read_cfg directories.results_dir)")"
PARALLEL_JOBS="$(read_cfg parameters.feature_parallel_jobs)"

IMPUTED_DIR="$RESULTS_DIR/Imputation/Output_Chromosomes/Beagle_1kG_imputed"
ORIGINAL_DIR="$RESULTS_DIR/Annotation/CADD_scored"
FEATURES_DIR="$RESULTS_DIR/Features"
FEATURE_SCRIPT="$REPO_DIR/Variant_to_gene/cal_features_multi.py"

if ! [[ "$PARALLEL_JOBS" =~ ^[0-9]+$ ]] || [[ "$PARALLEL_JOBS" -lt 1 ]]; then
	echo "ERROR: parameters.feature_parallel_jobs must be a positive integer." >&2
	exit 1
fi

if [[ ! -d "$ORIGINAL_DIR" ]]; then
	echo "ERROR: CADD-scored VCF directory not found: $ORIGINAL_DIR" >&2
	exit 1
fi

if [[ ! -d "$IMPUTED_DIR" ]]; then
	echo "WARNING: Imputed VCF directory not found: $IMPUTED_DIR" >&2
	echo "WARNING: Features will be calculated from original CADD-scored VCFs only." >&2
fi

if [[ ! -f "$FEATURE_SCRIPT" ]]; then
	echo "ERROR: Feature calculation script not found: $FEATURE_SCRIPT" >&2
	exit 1
fi

mkdir -p "$FEATURES_DIR"

mapfile -t CHROMOSOMES < <(
	find "$ORIGINAL_DIR" -maxdepth 1 -name "chr*.hg38_multianno.filtered_CADD_scored.vcf.gz" -printf "%f\n" \
	| sed -E 's/^chr//; s/\.hg38_multianno\.filtered_CADD_scored\.vcf\.gz$//' \
	| awk '$1 ~ /^([1-9]|1[0-9]|2[0-2]|X)$/ {print $1}' \
	| sort -V
)

if [[ "${#CHROMOSOMES[@]}" -eq 0 ]]; then
	echo "ERROR: No CADD-scored chromosome VCF files found in $ORIGINAL_DIR" >&2
	exit 1
fi

echo "Chromosomes selected for feature calculation: ${CHROMOSOMES[*]}"

calculate_features() {
	local chrom="$1"
	local original_vcf="$ORIGINAL_DIR/chr${chrom}.hg38_multianno.filtered_CADD_scored.vcf.gz"
	local imputed_vcf="$IMPUTED_DIR/chr${chrom}.filtered.vcf.gz"
	local out_dir="$FEATURES_DIR/chr${chrom}"

	if [[ -z "$chrom" ]]; then
		echo "ERROR: Empty chromosome argument in calculate_features" >&2
		return 1
	fi

	if [[ ! -f "$original_vcf" ]]; then
		echo "ERROR: Original CADD-scored VCF not found: $original_vcf" >&2
		return 1
	fi

	echo "Calculating features for chr${chrom}"
	mkdir -p "$out_dir"

	if [[ -f "$imputed_vcf" ]]; then
		python "$FEATURE_SCRIPT" \
			--original "$original_vcf" \
			--imputed "$imputed_vcf" \
			--out_dir "$out_dir" \
			2> >(grep -v '^\[W::bcf_hrec_check\] Invalid tag name:' >&2)
	else
		echo "WARNING: Filtered imputed VCF not found for chr${chrom}; calculating features from original VCF only." >&2

		python "$FEATURE_SCRIPT" \
			--original "$original_vcf" \
			--out_dir "$out_dir" \
			2> >(grep -v '^\[W::bcf_hrec_check\] Invalid tag name:' >&2)
	fi
}

run_jobs() {
	local fn="$1"
	shift

	local active_jobs=0
	local failed=0

	for chrom in "$@"; do
		"$fn" "$chrom" &
		active_jobs=$((active_jobs + 1))

		if [[ "$active_jobs" -ge "$PARALLEL_JOBS" ]]; then
			if ! wait -n; then
				failed=1
			fi
			active_jobs=$((active_jobs - 1))
		fi
	done

	while [[ "$active_jobs" -gt 0 ]]; do
		if ! wait -n; then
			failed=1
		fi
		active_jobs=$((active_jobs - 1))
	done

	if [[ "$failed" -ne 0 ]]; then
		echo "ERROR: At least one job failed during $fn" >&2
		exit 1
	fi
}

run_jobs calculate_features "${CHROMOSOMES[@]}"

echo "Gene-level feature calculation completed."