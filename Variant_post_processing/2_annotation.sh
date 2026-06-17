#!/usr/bin/env bash
set -euo pipefail

# This script performs post-imputation VCF processing:
# (1) annotates original chromosome VCFs with ANNOVAR,
# (2) retains functional variants,
# (3) intersects imputed VCFs with functionally filtered original variants,
# (4) annotates functionally filtered original VCFs with CADD scores.
# Usage: bash Variant_post_processing/2_annotation.sh config/local_config.yaml

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
ANNOVAR_DIR="$(python "$READ_CONFIG" "$CONFIG" tools.annovar_dir)"
CADD_TSV="$(python "$READ_CONFIG" "$CONFIG" resources.cadd_prescored)"
PARALLEL_JOBS="$(python "$READ_CONFIG" "$CONFIG" parameters.annotation_parallel_jobs)"

RESULTS_DIR="$(resolve_path "$RESULTS_DIR")"
ANNOVAR_DIR="$(resolve_path "$ANNOVAR_DIR")"
CADD_TSV="$(resolve_path "$CADD_TSV")"

TABLE_ANNOVAR="$ANNOVAR_DIR/table_annovar.pl"
ANNOVAR_DB="$ANNOVAR_DIR/humandb"

INPUT_DIR="$RESULTS_DIR/Imputation/Input_Chromosomes"
IMPUTED_DIR="$RESULTS_DIR/Imputation/Output_Chromosomes/Beagle_1kG_imputed"
ANNOTATION_DIR="$RESULTS_DIR/Annotation"
INTERMEDIATE_DIR="$RESULTS_DIR/Intermediate"
MULTIANNO_DIR="$ANNOTATION_DIR/multianno"
CADD_OUT_DIR="$ANNOTATION_DIR/CADD_scored"
HEADER="$ANNOTATION_DIR/cadd.header.txt"

ANNOVAR_BUILDVER="hg38"
ANNOVAR_PROTOCOLS="refGene,exac03,gnomad41_exome,ALL.sites.2015_08,avsnp150,clinvar_20241215,dbnsfp47a"
ANNOVAR_OPERATIONS="g,f,f,f,f,f,f"
BCFTOOLS_THREADS=3

CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X)

if [[ ! -d "$INPUT_DIR" ]]; then
	echo "ERROR: Input chromosome VCF directory not found: $INPUT_DIR"
	exit 1
fi

if [[ ! -d "$IMPUTED_DIR" ]]; then
	echo "ERROR: Imputed VCF directory not found: $IMPUTED_DIR"
	exit 1
fi

if [[ ! -f "$TABLE_ANNOVAR" ]]; then
	echo "ERROR: table_annovar.pl not found: $TABLE_ANNOVAR"
	exit 1
fi

if [[ ! -d "$ANNOVAR_DB" ]]; then
	echo "ERROR: ANNOVAR humandb directory not found: $ANNOVAR_DB"
	exit 1
fi

if [[ ! -f "$CADD_TSV" ]]; then
	echo "ERROR: CADD prescored file not found: $CADD_TSV"
	exit 1
fi

if [[ ! -f "${CADD_TSV}.tbi" ]]; then
	echo "ERROR: CADD tabix index not found: ${CADD_TSV}.tbi"
	exit 1
fi

mkdir -p "$MULTIANNO_DIR" "$CADD_OUT_DIR" "$INTERMEDIATE_DIR"

cat > "$HEADER" <<'HEADER_EOF'
##INFO=<ID=CADD1.7_RAW,Number=1,Type=Float,Description="CADD RawScore">
##INFO=<ID=CADD1.7_PHRED,Number=1,Type=Float,Description="CADD PHRED">
HEADER_EOF

annotate_with_annovar() {
	local chrom="$1"
	local in_vcf="$INPUT_DIR/chr${chrom}.vcf.gz"
	local out_prefix="$MULTIANNO_DIR/chr${chrom}"
	local multianno_vcf="${out_prefix}.${ANNOVAR_BUILDVER}_multianno.vcf"
	local filtered_vcf="$MULTIANNO_DIR/chr${chrom}.${ANNOVAR_BUILDVER}_multianno.function_filtered.vcf"

	echo "Annotating chr${chrom} with ANNOVAR"

	if [[ ! -f "$in_vcf" ]]; then
		echo "ERROR: Input VCF not found: $in_vcf"
		exit 1
	fi

	perl "$TABLE_ANNOVAR" \
		"$in_vcf" \
		"$ANNOVAR_DB" \
		-buildver "$ANNOVAR_BUILDVER" \
		-out "$out_prefix" \
		-remove \
		-protocol "$ANNOVAR_PROTOCOLS" \
		-operation "$ANNOVAR_OPERATIONS" \
		-nastring . \
		-polish \
		-vcfinput \
		-thread "$BCFTOOLS_THREADS"

	grep "^#" "$multianno_vcf" > "$filtered_vcf"

	grep -v "^#" "$multianno_vcf" | \
		grep -E "Func.refGene=(exonic|splicing|exonic;splicing|ncRNA_exonic;splicing)" \
		>> "$filtered_vcf" || true

	bgzip -f "$filtered_vcf"
	tabix -f -p vcf "${filtered_vcf}.gz"
}

filter_imputed_vcf() {
	local chrom="$1"
	local isec_dir="$INTERMEDIATE_DIR/chr${chrom}_isec_tmp"
	local imputed_vcf="$IMPUTED_DIR/chr${chrom}.vcf.gz"
	local filtered_original_vcf="$MULTIANNO_DIR/chr${chrom}.${ANNOVAR_BUILDVER}_multianno.function_filtered.vcf.gz"
	local out_vcf="$IMPUTED_DIR/chr${chrom}.filtered.vcf.gz"

	echo "Filtering imputed chr${chrom} to functional original variants"

	if [[ ! -f "$imputed_vcf" ]]; then
		echo "ERROR: Imputed VCF not found: $imputed_vcf"
		exit 1
	fi

	if [[ ! -f "$filtered_original_vcf" ]]; then
		echo "ERROR: Functionally filtered VCF not found: $filtered_original_vcf"
		exit 1
	fi

	rm -rf "$isec_dir"
	mkdir -p "$isec_dir"

	bcftools isec \
		-Oz \
		-p "$isec_dir" \
		-n=2 \
		-w1 \
		"$imputed_vcf" \
		"$filtered_original_vcf"

	mv "$isec_dir/0000.vcf.gz" "$out_vcf"
	bcftools index -f "$out_vcf"
}

annotate_with_cadd() {
	local chrom="$1"
	local in_vcf="$MULTIANNO_DIR/chr${chrom}.${ANNOVAR_BUILDVER}_multianno.function_filtered.vcf.gz"
	local out_vcf="$CADD_OUT_DIR/chr${chrom}.${ANNOVAR_BUILDVER}_multianno.filtered_CADD_scored.vcf.gz"

	echo "Annotating chr${chrom} with CADD"

	if [[ ! -f "$in_vcf" ]]; then
		echo "ERROR: Functionally filtered VCF not found: $in_vcf"
		exit 1
	fi

	bcftools annotate \
		-a "$CADD_TSV" \
		--threads "$BCFTOOLS_THREADS" \
		-c CHROM,POS,REF,ALT,INFO/CADD1.7_RAW:=5,INFO/CADD1.7_PHRED:=6 \
		-h "$HEADER" \
		-Oz \
		-o "$out_vcf" \
		"$in_vcf"

	bcftools index -f "$out_vcf"
}

export -f annotate_with_annovar filter_imputed_vcf annotate_with_cadd
export INPUT_DIR IMPUTED_DIR MULTIANNO_DIR CADD_OUT_DIR INTERMEDIATE_DIR
export TABLE_ANNOVAR ANNOVAR_DB CADD_TSV HEADER
export ANNOVAR_BUILDVER ANNOVAR_PROTOCOLS ANNOVAR_OPERATIONS BCFTOOLS_THREADS

printf "%s\n" "${CHROMOSOMES[@]}" | \
	parallel --halt soon,fail=1 -j "$PARALLEL_JOBS" \
	bash -c 'annotate_with_annovar "$1"' _ {}

echo "ANNOVAR annotation completed."

printf "%s\n" "${CHROMOSOMES[@]}" | \
	parallel --halt soon,fail=1 -j "$PARALLEL_JOBS" \
	bash -c 'filter_imputed_vcf "$1"' _ {}

echo "Imputed VCF filtering completed."

printf "%s\n" "${CHROMOSOMES[@]}" | \
	parallel --halt soon,fail=1 -j "$PARALLEL_JOBS" \
	bash -c 'annotate_with_cadd "$1"' _ {}

echo "CADD annotation completed."
