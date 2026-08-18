#!/usr/bin/env bash
set -euo pipefail

# This script performs genotype preprocessing and imputation:
# (1) Normalizes and filters a multisample VCF (split multiallelics, left-align, keep SNVs with AD > 5),
#     and splits data into per-chromosome VCF files
# (2) Conforms genotypes to a reference panel using conform-gt
# (3) Performs genotype imputation with Beagle using a population reference panel and genetic maps
#
# Input: multisample VCF (e.g. GLnexus output)
# Output: per-chromosome imputed and indexed VCF files

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

RESULTS_DIR="$(resolve_path "$(read_cfg directories.results_dir)")"
TMPDIR="$(resolve_path "$(read_cfg directories.temporary_dir)")"

REFERENCE_FASTA="$(resolve_path "$(read_cfg resources.reference_fasta)")"
REFERENCE_PANEL_DIR="$(resolve_path "$(read_cfg resources.reference_panel_dir)")"
MAP_DIR="$(resolve_path "$(read_cfg resources.genetic_maps_dir)")"

BEAGLE_JAR="$(resolve_path "$(read_cfg tools.beagle_jar)")"
CONFORM_JAR="$(resolve_path "$(read_cfg tools.conform_gt_jar)")"

PARALLEL_JOBS="$(read_cfg parameters.imputation_parallel_jobs)"
AD_THRESHOLD="$(read_cfg parameters.imputation_ad_threshold)"
BEAGLE_THREADS="$(read_cfg parameters.beagle_threads)"

GLNEXUS_OUTPUT="$RESULTS_DIR/Genotyping"
IMPUTATION_DIR="$RESULTS_DIR/Imputation"

IN_VCF="$GLNEXUS_OUTPUT/genotyped_multisample.vcf.gz"
INPUT_CHR_DIR="$IMPUTATION_DIR/Input_Chromosomes"
OUTPUT_CHR_DIR="$IMPUTATION_DIR/Output_Chromosomes/Beagle_1kG_imputed"

BCFTOOLS_THREADS=3
BEAGLE_MEMORY="50g"

if [[ ! -f "$IN_VCF" ]]; then
	echo "ERROR: Input multisample VCF not found: $IN_VCF" >&2
	exit 1
fi

if [[ ! -f "$REFERENCE_FASTA" ]]; then
	echo "ERROR: Reference FASTA not found: $REFERENCE_FASTA" >&2
	exit 1
fi

if [[ ! -f "$BEAGLE_JAR" ]]; then
	echo "ERROR: Beagle JAR not found: $BEAGLE_JAR" >&2
	exit 1
fi

if [[ ! -f "$CONFORM_JAR" ]]; then
	echo "ERROR: conform-gt JAR not found: $CONFORM_JAR" >&2
	exit 1
fi

if [[ ! -d "$MAP_DIR" ]]; then
	echo "ERROR: Genetic maps directory not found: $MAP_DIR" >&2
	exit 1
fi

if [[ ! -d "$REFERENCE_PANEL_DIR" ]]; then
	echo "ERROR: Reference panel directory not found: $REFERENCE_PANEL_DIR" >&2
	exit 1
fi

mkdir -p "$INPUT_CHR_DIR" "$OUTPUT_CHR_DIR" "$TMPDIR"

if [[ ! -f "${IN_VCF}.tbi" && ! -f "${IN_VCF}.csi" ]]; then
	echo "Indexing input multisample VCF"
	bcftools index --threads "$BCFTOOLS_THREADS" -f "$IN_VCF"
fi

mapfile -t CHROMOSOMES < <(
	bcftools index -s "$IN_VCF" \
	| cut -f1 \
	| awk '$1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1}' \
	| sort -V
)

if [[ "${#CHROMOSOMES[@]}" -eq 0 ]]; then
	echo "ERROR: No canonical chromosomes found in $IN_VCF" >&2
	exit 1
fi

echo "Chromosomes selected from input VCF: ${CHROMOSOMES[*]}"

preprocess_chromosome() {
	local chrom="$1"
	local out_vcf="$INPUT_CHR_DIR/chr${chrom}.vcf.gz"

	if [[ -z "$chrom" ]]; then
		echo "ERROR: Empty chromosome argument in preprocess_chromosome" >&2
		return 1
	fi

	if [[ -s "$out_vcf" && -s "${out_vcf}.tbi" ]]; then
		echo "Skipping preprocessing chr${chrom} - output exists"
		return 0
	fi

	echo "Preprocessing chr${chrom}"

	bcftools norm \
		-m -both \
		-f "$REFERENCE_FASTA" \
		-d all \
		-r "$chrom" \
		--threads "$BCFTOOLS_THREADS" \
		-Oz "$IN_VCF" \
	| bcftools view \
		-v snps \
		-i "MAX(FMT/AD[*]) > ${AD_THRESHOLD}" \
		--threads "$BCFTOOLS_THREADS" \
		-Oz \
		-o "$out_vcf"

	bcftools index \
		--threads "$BCFTOOLS_THREADS" \
		-t -f "$out_vcf"
}

conform_chromosome() {
	local chrom="$1"
	local gt_vcf="$INPUT_CHR_DIR/chr${chrom}.vcf.gz"
	local ref_vcf
	ref_vcf=$(find "$REFERENCE_PANEL_DIR" -maxdepth 1 -type f -name "1kGP.${chrom}.*.nochr.vcf.gz" | head -n 1)
	local out_prefix="$INPUT_CHR_DIR/chr${chrom}.conformed"
	local out_vcf="${out_prefix}.vcf.gz"

	if [[ -z "$chrom" ]]; then
		echo "ERROR: Empty chromosome argument in conform_chromosome" >&2
		return 1
	fi

	if [[ -s "$out_vcf" && -s "${out_vcf}.tbi" ]]; then
		echo "Skipping conformation chr${chrom} - output exists"
		return 0
	fi

	echo "Conforming chr${chrom}"

	if [[ ! -f "$gt_vcf" ]]; then
		echo "ERROR: Input chromosome VCF not found: $gt_vcf" >&2
		return 1
	fi

	if [[ ! -f "$ref_vcf" ]]; then
		echo "ERROR: Reference panel VCF not found for chr${chrom}: $ref_vcf" >&2
		return 1
	fi

	java -jar "$CONFORM_JAR" \
		gt="$gt_vcf" \
		ref="$ref_vcf" \
		chrom="$chrom" \
		match=POS \
		out="$out_prefix"

	bcftools index \
		--threads "$BCFTOOLS_THREADS" \
		-t -f "$out_vcf"
}

impute_chromosome() {
	local chrom="$1"
	local gt_vcf="$INPUT_CHR_DIR/chr${chrom}.conformed.vcf.gz"
	local map_file="$MAP_DIR/plink.chr${chrom}.GRCh38.map"
	local ref_bref
	ref_bref=$(find "$REFERENCE_PANEL_DIR/BREF3" -maxdepth 1 -type f -name "1kGP.${chrom}.*.nochr.bref3" | head -n 1)
	local out_prefix="$OUTPUT_CHR_DIR/chr${chrom}"
	local out_vcf="${out_prefix}.vcf.gz"

	if [[ -z "$chrom" ]]; then
		echo "ERROR: Empty chromosome argument in impute_chromosome" >&2
		return 1
	fi

	if [[ -s "$out_vcf" && -s "${out_vcf}.tbi" ]]; then
		echo "Skipping imputation chr${chrom} - output exists"
		return 0
	fi

	echo "Imputing chr${chrom}"

	if [[ ! -f "$gt_vcf" ]]; then
		echo "ERROR: Conformed VCF not found: $gt_vcf" >&2
		return 1
	fi

	if [[ ! -f "$map_file" ]]; then
		echo "ERROR: Genetic map file not found for chr${chrom}: $map_file" >&2
		return 1
	fi

	if [[ ! -f "$ref_bref" ]]; then
		echo "ERROR: Reference panel BREF3 file not found for chr${chrom}: $ref_bref" >&2
		return 1
	fi

	java -Xmx"$BEAGLE_MEMORY" -jar "$BEAGLE_JAR" \
		gt="$gt_vcf" \
		map="$map_file" \
		ref="$ref_bref" \
		nthreads="$BEAGLE_THREADS" \
		out="$out_prefix"

	bcftools index \
		--threads "$BCFTOOLS_THREADS" \
		-t -f "$out_vcf"
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

run_jobs preprocess_chromosome "${CHROMOSOMES[@]}"

echo "Preprocessing completed."

mapfile -t CHROMOSOMES_AFTER_FILTERING < <(
	for chrom in "${CHROMOSOMES[@]}"; do
		vcf="$INPUT_CHR_DIR/chr${chrom}.vcf.gz"

		if [[ -f "$vcf" && -f "${vcf}.tbi" ]]; then
			records="$(bcftools index -n "$vcf" 2>/dev/null || echo 0)"

			if [[ "$records" -gt 0 ]]; then
				echo "$chrom"
			else
				echo "WARNING: No SNV records retained for chr${chrom} after filtering; skipping downstream imputation." >&2
			fi
		fi
	done
)

if [[ "${#CHROMOSOMES_AFTER_FILTERING[@]}" -eq 0 ]]; then
	echo "ERROR: No chromosomes retained after preprocessing and AD filtering." >&2
	exit 1
fi

echo "Chromosomes retained after filtering: ${CHROMOSOMES_AFTER_FILTERING[*]}"

run_jobs conform_chromosome "${CHROMOSOMES_AFTER_FILTERING[@]}"

echo "Genotype conformation completed."

run_jobs impute_chromosome "${CHROMOSOMES_AFTER_FILTERING[@]}"

echo "Genotype imputation completed."
