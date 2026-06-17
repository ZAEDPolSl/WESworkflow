#!/usr/bin/env bash
set -euo pipefail

# This script runs GLnexus to merge and jointly genotype GVCF files
# produced by DeepVariant, then converts the output BCF to VCF.
# Usage: bash Variant_calling/Joint_genotyping/run_GLnexus.sh config/local_config.yaml

CONFIG="${1:-config/local_config.yaml}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

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

GVCF_SEARCH_ROOT="$(python "$READ_CONFIG" "$CONFIG" directories.gvcf_search_root)"
RESULTS_DIR="$(python "$READ_CONFIG" "$CONFIG" directories.results_dir)"
RESULTS_DIR="$(resolve_path "$RESULTS_DIR")"
OUTPUT_DIR="$RESULTS_DIR/Genotyping"
TMPDIR="$(python "$READ_CONFIG" "$CONFIG" directories.temporary_dir)"
REGIONS_BED="$(python "$READ_CONFIG" "$CONFIG" resources.refseq_bed)"
BCFTOOLS_THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.bcftools_threads)"

REGIONS_BED="$(resolve_path "$REGIONS_BED")"
SCRATCH_DIR="$TMPDIR/GLnexus_scratch"
GLNEXUS_CONFIG="DeepVariantWES"

if [[ ! -d "$GVCF_SEARCH_ROOT" ]]; then
	echo "ERROR: GVCF search root not found: $GVCF_SEARCH_ROOT"
	exit 1
fi

if [[ ! -f "$REGIONS_BED" ]]; then
	echo "ERROR: Regions BED file not found: $REGIONS_BED"
	exit 1
fi

mkdir -p "$OUTPUT_DIR" "$TMPDIR"
rm -rf "$SCRATCH_DIR"

readarray -d '' GVCF_FILES < <(
	find "$GVCF_SEARCH_ROOT" -type f -name "*.g.vcf.gz" -print0 | sort -z
)

if [[ "${#GVCF_FILES[@]}" -eq 0 ]]; then
	echo "ERROR: No *.g.vcf.gz files found under $GVCF_SEARCH_ROOT"
	exit 1
fi

echo "Found ${#GVCF_FILES[@]} GVCF files."
echo "Running GLnexus joint genotyping."

glnexus_cli \
	--config "$GLNEXUS_CONFIG" \
	--bed "$REGIONS_BED" \
	--dir "$SCRATCH_DIR" \
	"${GVCF_FILES[@]}" \
	> "$OUTPUT_DIR/genotyped_multisample.bcf"

echo "Converting BCF to uncompressed VCF."

bcftools view \
	--threads "$BCFTOOLS_THREADS" \
	-Ov "$OUTPUT_DIR/genotyped_multisample.bcf" \
	-o "$OUTPUT_DIR/genotyped_multisample.vcf"

echo "Converting BCF to compressed VCF."

bcftools view \
	--threads "$BCFTOOLS_THREADS" \
	-Oz "$OUTPUT_DIR/genotyped_multisample.bcf" \
	-o "$OUTPUT_DIR/genotyped_multisample.vcf.gz"

echo "Indexing compressed VCF."

bcftools index \
	--threads "$BCFTOOLS_THREADS" \
	-f "$OUTPUT_DIR/genotyped_multisample.vcf.gz"

echo "GLnexus joint genotyping completed successfully."
