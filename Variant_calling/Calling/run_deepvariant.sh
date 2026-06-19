#!/usr/bin/env bash
set -euo pipefail

# This script performs variant calling on WES BAM files using DeepVariant
# in a Docker container, producing per-sample VCF and GVCF files.
# Usage: bash Variant_calling/Calling/run_deepvariant.sh config/local_config.yaml

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
	echo "ERROR: Config file not found: $CONFIG" >&2
	exit 1
fi

BAM_DIR="$(resolve_path "$(python "$READ_CONFIG" "$CONFIG" directories.deepvariant_bam)")"
OUTPUT_DIR="$(resolve_path "$(python "$READ_CONFIG" "$CONFIG" directories.deepvariant_output)")"
TMPDIR="$(resolve_path "$(python "$READ_CONFIG" "$CONFIG" directories.temporary_dir)")"

REFERENCE_FASTA="$(resolve_path "$(python "$READ_CONFIG" "$CONFIG" resources.reference_fasta)")"
REGIONS_BED="$(resolve_path "$(python "$READ_CONFIG" "$CONFIG" resources.refseq_bed)")"

THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.deepvariant_threads)"
USE_GPU="$(python "$READ_CONFIG" "$CONFIG" parameters.deepvariant_use_gpu | tr '[:upper:]' '[:lower:]')"
BAM_SUFFIX="$(python "$READ_CONFIG" "$CONFIG" directories.deepvariant_bam_suffix 2>/dev/null || echo "_marked.bam")"
DOCKER_IMAGE="$(python "$READ_CONFIG" "$CONFIG" tools.deepvariant_docker_image)"

MODEL_TYPE="WES"

REFERENCE_DIR="$(dirname "$REFERENCE_FASTA")"
REFERENCE_FILE="$(basename "$REFERENCE_FASTA")"
REGIONS_DIR="$(dirname "$REGIONS_BED")"
REGIONS_FILE="$(basename "$REGIONS_BED")"

if [[ ! -d "$BAM_DIR" ]]; then
	echo "ERROR: Input BAM directory not found: $BAM_DIR" >&2
	exit 1
fi

if [[ ! -f "$REFERENCE_FASTA" ]]; then
	echo "ERROR: Reference FASTA not found: $REFERENCE_FASTA" >&2
	exit 1
fi

if [[ ! -f "$REGIONS_BED" ]]; then
	echo "ERROR: Regions BED file not found: $REGIONS_BED" >&2
	exit 1
fi

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/logs" "$TMPDIR"

DOCKER_GPU_ARGS=()
if [[ "$USE_GPU" == "true" ]]; then
	DOCKER_GPU_ARGS=(--gpus 1)
fi

total="$(find "$BAM_DIR" -name "*${BAM_SUFFIX}" | wc -l)"

if [[ "$total" -eq 0 ]]; then
	echo "ERROR: No BAM files found in $BAM_DIR using suffix $BAM_SUFFIX" >&2
	exit 1
fi

i=1
failed_samples=()

while IFS= read -r -d '' bam_file; do
	filename="$(basename "$bam_file")"
	sample_name="${filename%"$BAM_SUFFIX"}"
	relative_bam="${bam_file#"$BAM_DIR"/}"

	echo "[$i/$total] Processing $sample_name"

	if [[ -f "$OUTPUT_DIR/${sample_name}.vcf.gz" && -f "$OUTPUT_DIR/${sample_name}.g.vcf.gz" ]]; then
		echo "Skipping $sample_name - VCF and GVCF exist"
		i=$((i + 1))
		continue
	fi

	if ! docker run --rm "${DOCKER_GPU_ARGS[@]}" \
		-v "$BAM_DIR":"/input" \
		-v "$OUTPUT_DIR":"/output" \
		-v "$REFERENCE_DIR":"/reference" \
		-v "$REGIONS_DIR":"/regions" \
		-v "$TMPDIR":"/tmp" \
		"$DOCKER_IMAGE" \
		/opt/deepvariant/bin/run_deepvariant \
		--model_type "$MODEL_TYPE" \
		--ref "/reference/$REFERENCE_FILE" \
		--reads "/input/$relative_bam" \
		--regions "/regions/$REGIONS_FILE" \
		--output_vcf "/output/${sample_name}.vcf.gz" \
		--output_gvcf "/output/${sample_name}.g.vcf.gz" \
		--logging_dir "/output/logs" \
		--num_shards "$THREADS" \
		--intermediate_results_dir "/tmp/dv_${sample_name}" \
		--sample_name "$sample_name"; then
		echo "ERROR: DeepVariant failed for sample $sample_name" >&2
		failed_samples+=("$sample_name")
	fi

	i=$((i + 1))
done < <(find "$BAM_DIR" -name "*${BAM_SUFFIX}" -print0 | sort -z)

if [[ "${#failed_samples[@]}" -gt 0 ]]; then
	echo "The following samples failed:" >&2
	printf '%s\n' "${failed_samples[@]}" >&2
	exit 1
fi

echo "DeepVariant completed successfully."