#!/usr/bin/env bash
set -euo pipefail

# This script performs genotype preprocessing and imputation:
# (1) normalizes and filters a multisample VCF,
# (2) splits data into per-chromosome VCF files,
# (3) conforms genotypes to a reference panel using conform-gt,
# (4) performs genotype imputation with Beagle.
# Usage: bash Variant_post_processing/1_genotype_imputation.sh config/local_config.yaml

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
TMPDIR="$(python "$READ_CONFIG" "$CONFIG" directories.temporary_dir)"

REFERENCE_FASTA="$(python "$READ_CONFIG" "$CONFIG" resources.reference_fasta)"
REFERENCE_PANEL_DIR="$(python "$READ_CONFIG" "$CONFIG" resources.reference_panel_dir)"
MAP_DIR="$(python "$READ_CONFIG" "$CONFIG" resources.genetic_maps_dir)"

BEAGLE_JAR="$(python "$READ_CONFIG" "$CONFIG" tools.beagle_jar)"
CONFORM_JAR="$(python "$READ_CONFIG" "$CONFIG" tools.conform_gt_jar)"

PARALLEL_JOBS="$(python "$READ_CONFIG" "$CONFIG" parameters.imputation_parallel_jobs)"
AD_THRESHOLD="$(python "$READ_CONFIG" "$CONFIG" parameters.imputation_ad_threshold)"
BEAGLE_THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.beagle_threads)"

RESULTS_DIR="$(resolve_path "$RESULTS_DIR")"
TMPDIR="$(resolve_path "$TMPDIR")"
REFERENCE_FASTA="$(resolve_path "$REFERENCE_FASTA")"
REFERENCE_PANEL_DIR="$(resolve_path "$REFERENCE_PANEL_DIR")"
MAP_DIR="$(resolve_path "$MAP_DIR")"
BEAGLE_JAR="$(resolve_path "$BEAGLE_JAR")"
CONFORM_JAR="$(resolve_path "$CONFORM_JAR")"

GLNEXUS_OUTPUT="$RESULTS_DIR/Genotyping"
IMPUTATION_DIR="$RESULTS_DIR/Imputation"

IN_VCF="$GLNEXUS_OUTPUT/genotyped_multisample.vcf.gz"
INPUT_CHR_DIR="$IMPUTATION_DIR/Input_Chromosomes"
OUTPUT_CHR_DIR="$IMPUTATION_DIR/Output_Chromosomes/Beagle_1kG_imputed"

BCFTOOLS_THREADS=3
BEAGLE_MEMORY="50g"

CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

if [[ ! -f "$IN_VCF" ]]; then
	echo "ERROR: Input multisample VCF not found: $IN_VCF"
	exit 1
fi

if [[ ! -f "$REFERENCE_FASTA" ]]; then
	echo "ERROR: Reference FASTA not found: $REFERENCE_FASTA"
	exit 1
fi

if [[ ! -f "$BEAGLE_JAR" ]]; then
	echo "ERROR: Beagle JAR not found: $BEAGLE_JAR"
	exit 1
fi

if [[ ! -f "$CONFORM_JAR" ]]; then
	echo "ERROR: conform-gt JAR not found: $CONFORM_JAR"
	exit 1
fi

if [[ ! -d "$MAP_DIR" ]]; then
	echo "ERROR: Genetic maps directory not found: $MAP_DIR"
	exit 1
fi

if [[ ! -d "$REFERENCE_PANEL_DIR" ]]; then
	echo "ERROR: Reference panel directory not found: $REFERENCE_PANEL_DIR"
	exit 1
fi

mkdir -p "$INPUT_CHR_DIR" "$OUTPUT_CHR_DIR" "$TMPDIR"

preprocess_chromosome() {
	local chrom="$1"
	local out_vcf="$INPUT_CHR_DIR/chr${chrom}.vcf.gz"

	echo "Preprocessing chr${chrom}"

	bcftools norm \
		-m -both \
		-f "$REFERENCE_FASTA" \
		-d all \
		-r "$chrom" \
		--threads "$BCFTOOLS_THREADS" \
		-Oz "$IN_VCF" | \
	bcftools view \
		-v snps \
		-i "MAX(FMT/AD[*]) > ${AD_THRESHOLD}" \
		--threads "$BCFTOOLS_THREADS" \
		-Oz \
		-o "$out_vcf"

	bcftools index \
		--threads "$BCFTOOLS_THREADS" \
		-f "$out_vcf"
}

conform_chromosome() {
	local chrom="$1"
	local gt_vcf="$INPUT_CHR_DIR/chr${chrom}.vcf.gz"
	local ref_vcf="$REFERENCE_PANEL_DIR/1kGP.${chrom}.EUR.nochr.vcf.gz"
	local out_prefix="$INPUT_CHR_DIR/chr${chrom}.conformed"
	local out_vcf="${out_prefix}.vcf.gz"

	echo "Conforming chr${chrom}"

	if [[ ! -f "$gt_vcf" ]]; then
		echo "ERROR: Input chromosome VCF not found: $gt_vcf"
		exit 1
	fi

	if [[ ! -f "$ref_vcf" ]]; then
		echo "ERROR: Reference panel VCF not found: $ref_vcf"
		exit 1
	fi

	java -jar "$CONFORM_JAR" \
		gt="$gt_vcf" \
		ref="$ref_vcf" \
		chrom="$chrom" \
		match=POS \
		out="$out_prefix"

	bcftools index \
		--threads "$BCFTOOLS_THREADS" \
		-f "$out_vcf"
}

impute_chromosome() {
	local chrom="$1"
	local gt_vcf="$INPUT_CHR_DIR/chr${chrom}.conformed.vcf.gz"
	local map_file="$MAP_DIR/plink.chr${chrom}.GRCh38.map"
	local ref_bref="$REFERENCE_PANEL_DIR/BREF3/1kGP.${chrom}.EUR.nochr.bref3"
	local out_prefix="$OUTPUT_CHR_DIR/chr${chrom}"
	local out_vcf="${out_prefix}.vcf.gz"

	echo "Imputing chr${chrom}"

	if [[ ! -f "$gt_vcf" ]]; then
		echo "ERROR: Conformed VCF not found: $gt_vcf"
		exit 1
	fi

	if [[ ! -f "$map_file" ]]; then
		echo "ERROR: Genetic map file not found: $map_file"
		exit 1
	fi

	if [[ ! -f "$ref_bref" ]]; then
		echo "ERROR: Reference panel BREF3 file not found: $ref_bref"
		exit 1
	fi

	java -Xmx"$BEAGLE_MEMORY" -jar "$BEAGLE_JAR" \
		gt="$gt_vcf" \
		map="$map_file" \
		ref="$ref_bref" \
		nthreads="$BEAGLE_THREADS" \
		out="$out_prefix"

	bcftools index \
		--threads "$BCFTOOLS_THREADS" \
		-f "$out_vcf"
}

export -f preprocess_chromosome conform_chromosome impute_chromosome
export IN_VCF INPUT_CHR_DIR OUTPUT_CHR_DIR
export REFERENCE_FASTA REFERENCE_PANEL_DIR MAP_DIR
export BEAGLE_JAR CONFORM_JAR
export BCFTOOLS_THREADS AD_THRESHOLD BEAGLE_THREADS BEAGLE_MEMORY

printf "%s\n" "${CHROMOSOMES[@]}" | \
	parallel --halt soon,fail=1 -j "$PARALLEL_JOBS" \
	bash -c 'preprocess_chromosome "$1"' _ {}

echo "Preprocessing completed."

printf "%s\n" "${CHROMOSOMES[@]}" | \
	parallel --line-buffer --halt soon,fail=1 -j "$PARALLEL_JOBS" \
	bash -c 'conform_chromosome "$1"' _ {}

echo "Genotype conformation completed."

printf "%s\n" "${CHROMOSOMES[@]}" | \
	parallel --halt soon,fail=1 -j "$PARALLEL_JOBS" \
	bash -c 'impute_chromosome "$1"' _ {}

echo "Genotype imputation completed."
