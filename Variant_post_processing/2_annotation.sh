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

read_cfg() {
	python "$READ_CONFIG" "$CONFIG" "$1"
}

if [[ ! -f "$CONFIG" ]]; then
	echo "ERROR: Config file not found: $CONFIG" >&2
	exit 1
fi

RESULTS_DIR="$(resolve_path "$(read_cfg directories.results_dir)")"
ANNOVAR_DIR="$(resolve_path "$(read_cfg tools.annovar_dir)")"
CADD_TSV="$(resolve_path "$(read_cfg resources.cadd_prescored)")"
PARALLEL_JOBS="$(read_cfg parameters.annotation_parallel_jobs)"

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

if [[ ! -d "$INPUT_DIR" ]]; then
	echo "ERROR: Input chromosome VCF directory not found: $INPUT_DIR" >&2
	exit 1
fi

if [[ ! -d "$IMPUTED_DIR" ]]; then
	echo "ERROR: Imputed VCF directory not found: $IMPUTED_DIR" >&2
	exit 1
fi

if [[ ! -f "$TABLE_ANNOVAR" ]]; then
	echo "ERROR: table_annovar.pl not found: $TABLE_ANNOVAR" >&2
	exit 1
fi

if [[ ! -d "$ANNOVAR_DB" ]]; then
	echo "ERROR: ANNOVAR humandb directory not found: $ANNOVAR_DB" >&2
	exit 1
fi

if [[ ! -f "$CADD_TSV" ]]; then
	echo "ERROR: CADD prescored file not found: $CADD_TSV" >&2
	exit 1
fi

if [[ ! -f "${CADD_TSV}.tbi" ]]; then
	echo "ERROR: CADD tabix index not found: ${CADD_TSV}.tbi" >&2
	exit 1
fi

mkdir -p "$MULTIANNO_DIR" "$CADD_OUT_DIR" "$INTERMEDIATE_DIR"

cat > "$HEADER" <<'HEADER_EOF'
##INFO=<ID=CADD1.7_RAW,Number=1,Type=Float,Description="CADD RawScore">
##INFO=<ID=CADD1.7_PHRED,Number=1,Type=Float,Description="CADD PHRED">
HEADER_EOF

vcf_has_records() {
	local vcf="$1"
	local n

	[[ -f "$vcf" ]] || return 1
	[[ -f "${vcf}.tbi" || -f "${vcf}.csi" ]] || return 1

	n="$(bcftools index -n "$vcf" 2>/dev/null || echo 0)"
	[[ "$n" =~ ^[0-9]+$ && "$n" -gt 0 ]]
}

get_input_chromosomes() {
	find "$INPUT_DIR" -maxdepth 1 -name "chr*.vcf.gz" -printf "%f\n" \
	| sed -E 's/^chr//; s/\.vcf\.gz$//' \
	| awk '$1 ~ /^([1-9]|1[0-9]|2[0-2]|X|Y)$/ {print $1}' \
	| sort -V
}

mapfile -t CHROMOSOMES < <(get_input_chromosomes)

if [[ "${#CHROMOSOMES[@]}" -eq 0 ]]; then
	echo "ERROR: No chromosome VCF files found in $INPUT_DIR" >&2
	exit 1
fi

echo "Chromosomes selected for annotation: ${CHROMOSOMES[*]}"

annotate_with_annovar() {
	local chrom="$1"
	local in_vcf="$INPUT_DIR/chr${chrom}.vcf.gz"
	local out_prefix="$MULTIANNO_DIR/chr${chrom}"
	local multianno_vcf="${out_prefix}.${ANNOVAR_BUILDVER}_multianno.vcf"
	local filtered_vcf="$MULTIANNO_DIR/chr${chrom}.${ANNOVAR_BUILDVER}_multianno.function_filtered.vcf"

	if [[ -z "$chrom" ]]; then
		echo "ERROR: Empty chromosome argument in annotate_with_annovar" >&2
		return 1
	fi

	if [[ -f "${filtered_vcf}.gz" && -f "${filtered_vcf}.gz.tbi" ]]; then
		echo "Skipping ANNOVAR chr${chrom} - output exists"
		return 0
	fi

	echo "Annotating chr${chrom} with ANNOVAR"

	if [[ ! -f "$in_vcf" ]]; then
		echo "ERROR: Input VCF not found: $in_vcf" >&2
		return 1
	fi

	if ! vcf_has_records "$in_vcf"; then
		echo "WARNING: Input VCF has no records for chr${chrom}; skipping ANNOVAR." >&2
		return 0
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

	grep -v "^#" "$multianno_vcf" \
	| grep -E "Func.refGene=(exonic|splicing|exonic;splicing|ncRNA_exonic;splicing)" \
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

	if [[ -z "$chrom" ]]; then
		echo "ERROR: Empty chromosome argument in filter_imputed_vcf" >&2
		return 1
	fi

	if [[ -f "$out_vcf" && -f "${out_vcf}.tbi" ]]; then
		echo "Skipping imputed filtering chr${chrom} - output exists"
		return 0
	fi

	echo "Filtering imputed chr${chrom} to functional original variants"

	if [[ ! -f "$imputed_vcf" ]]; then
		echo "WARNING: Imputed VCF not found for chr${chrom}; skipping imputed filtering." >&2
		return 0
	fi

	if [[ ! -f "$filtered_original_vcf" ]]; then
		echo "WARNING: Functionally filtered VCF not found for chr${chrom}; skipping imputed filtering." >&2
		return 0
	fi

	if ! vcf_has_records "$filtered_original_vcf"; then
		echo "WARNING: No functional original variants for chr${chrom}; skipping imputed filtering." >&2
		return 0
	fi

	if ! vcf_has_records "$imputed_vcf"; then
		echo "WARNING: Imputed VCF has no records for chr${chrom}; skipping imputed filtering." >&2
		return 0
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

	if [[ -f "$isec_dir/0000.vcf.gz" ]]; then
		mv "$isec_dir/0000.vcf.gz" "$out_vcf"
		bcftools index -f "$out_vcf"
	else
		echo "WARNING: No intersection output for chr${chrom}; skipping." >&2
	fi
}

annotate_with_cadd() {
	local chrom="$1"
	local in_vcf="$MULTIANNO_DIR/chr${chrom}.${ANNOVAR_BUILDVER}_multianno.function_filtered.vcf.gz"
	local out_vcf="$CADD_OUT_DIR/chr${chrom}.${ANNOVAR_BUILDVER}_multianno.filtered_CADD_scored.vcf.gz"

	if [[ -z "$chrom" ]]; then
		echo "ERROR: Empty chromosome argument in annotate_with_cadd" >&2
		return 1
	fi

	if [[ -f "$out_vcf" && -f "${out_vcf}.tbi" ]]; then
		echo "Skipping CADD chr${chrom} - output exists"
		return 0
	fi

	echo "Annotating chr${chrom} with CADD"

	if [[ ! -f "$in_vcf" ]]; then
		echo "WARNING: Functionally filtered VCF not found for chr${chrom}; skipping CADD." >&2
		return 0
	fi

	if ! vcf_has_records "$in_vcf"; then
		echo "WARNING: No functional variants for chr${chrom}; skipping CADD." >&2
		return 0
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

run_jobs annotate_with_annovar "${CHROMOSOMES[@]}"

echo "ANNOVAR annotation completed."

mapfile -t CHROMOSOMES_WITH_IMPUTED < <(
	for chrom in "${CHROMOSOMES[@]}"; do
		if [[ -f "$IMPUTED_DIR/chr${chrom}.vcf.gz" ]]; then
			echo "$chrom"
		else
			echo "WARNING: No imputed VCF for chr${chrom}; it will be skipped during imputed VCF filtering." >&2
		fi
	done
)

if [[ "${#CHROMOSOMES_WITH_IMPUTED[@]}" -gt 0 ]]; then
	echo "Chromosomes selected for imputed VCF filtering: ${CHROMOSOMES_WITH_IMPUTED[*]}"
	run_jobs filter_imputed_vcf "${CHROMOSOMES_WITH_IMPUTED[@]}"
else
	echo "WARNING: No imputed VCF files found. Skipping imputed VCF filtering." >&2
fi

echo "Imputed VCF filtering completed."

run_jobs annotate_with_cadd "${CHROMOSOMES[@]}"

echo "CADD annotation completed."