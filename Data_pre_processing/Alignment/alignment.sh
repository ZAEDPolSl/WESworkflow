#!/usr/bin/env bash
set -euo pipefail

# This script aligns paired-end FASTQ files to a reference genome using BWA MEM,
# sorts the resulting BAM files, and marks sequencing duplicates with Picard.
# Usage: bash Data_pre_processing/Alignment/alignment.sh config/local_config.yaml

CONFIG="${1:-config/local_config.yaml}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if [[ "$CONFIG" != /* ]]; then
	CONFIG="$REPO_DIR/$CONFIG"
fi

READ_CONFIG="$REPO_DIR/scripts/read_config.py"

if [[ ! -f "$CONFIG" ]]; then
	echo "ERROR: Config file not found: $CONFIG"
	exit 1
fi

THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.alignment_threads)"
BWA_MIN_SCORE="$(python "$READ_CONFIG" "$CONFIG" parameters.bwa_min_score)"
MAX_JOBS="$(python "$READ_CONFIG" "$CONFIG" parameters.mark_duplicates_jobs)"
R1_SUFFIX="$(python "$READ_CONFIG" "$CONFIG" parameters.fastq_r1_suffix)"
R2_SUFFIX="$(python "$READ_CONFIG" "$CONFIG" parameters.fastq_r2_suffix)"

REFERENCE_FASTA="$(python "$READ_CONFIG" "$CONFIG" resources.reference_fasta)"
INPUT_DIR="$(python "$READ_CONFIG" "$CONFIG" directories.alignment_fastq)"
OUTPUT_DIR="$(python "$READ_CONFIG" "$CONFIG" directories.bam)"
TMPDIR="$(python "$READ_CONFIG" "$CONFIG" directories.temporary_dir)"

VALIDATION_STRINGENCY="STRICT"
RG_LIBRARY="WES"
RG_PLATFORM="Illumina"

export TMPDIR

if [[ ! -f "$REFERENCE_FASTA" ]]; then
	echo "ERROR: Reference FASTA not found: $REFERENCE_FASTA"
	exit 1
fi

if [[ ! -d "$INPUT_DIR" ]]; then
	echo "ERROR: Input FASTQ directory not found: $INPUT_DIR"
	exit 1
fi

mkdir -p "$OUTPUT_DIR" "$TMPDIR"

echo "STARTED ALIGNMENT"

failed_samples=()
i=1
total="$(find "$INPUT_DIR" -name "*${R1_SUFFIX}" | wc -l)"

if [[ "$total" -eq 0 ]]; then
	echo "ERROR: No R1 FASTQ files found in $INPUT_DIR using suffix $R1_SUFFIX"
	exit 1
fi

while IFS= read -r -d '' r1; do
	filename="$(basename "$r1")"
	sample="${filename%"$R1_SUFFIX"}"
	r2="$INPUT_DIR/${sample}${R2_SUFFIX}"

	if [[ ! -f "$r2" ]]; then
		echo "WARNING: Missing R2 file for sample $sample: $r2"
		failed_samples+=("$sample")
		continue
	fi

	if [[ -f "$OUTPUT_DIR/${sample}.bam" ]]; then
		echo "Skipping $sample - BAM exists"
		((i++))
		continue
	fi

	echo "Processing sample $sample ($i/$total)"

	if ! bwa mem -t "$THREADS" -T "$BWA_MIN_SCORE" -v 1 \
		-R "@RG\tID:${sample}\tSM:${sample}\tLB:${RG_LIBRARY}\tPL:${RG_PLATFORM}" \
		"$REFERENCE_FASTA" "$r1" "$r2" | \
		samtools view -b -@ "$THREADS" - | \
		samtools sort -@ "$THREADS" -T "$TMPDIR/${sample}_tmp" \
			-o "$OUTPUT_DIR/${sample}.bam"; then

		echo "ERROR: Alignment or sorting failed for sample $sample"
		failed_samples+=("$sample")
	fi

	((i++))
done < <(find "$INPUT_DIR" -name "*${R1_SUFFIX}" -print0 | sort -z)

if [[ "${#failed_samples[@]}" -gt 0 ]]; then
	echo "The following samples failed during alignment:"
	printf '%s\n' "${failed_samples[@]}"
else
	echo "ALIGNMENT COMPLETED SUCCESSFULLY"
fi

echo "STARTED DUPLICATE MARKING"

mark_duplicates() {
	local bam_file="$1"
	local output_bam
	local metrics_file

	output_bam="${bam_file%.bam}_marked.bam"
	metrics_file="${bam_file%.bam}_metrics.txt"

	if [[ -e "$output_bam" || "$bam_file" == *_marked.bam ]]; then
		echo "Skipping $bam_file"
		return
	fi

	echo "Marking duplicates in $bam_file"

	if ! picard MarkDuplicates \
		--INPUT "$bam_file" \
		--OUTPUT "$output_bam" \
		--METRICS_FILE "$metrics_file" \
		--CREATE_INDEX true \
		--VALIDATION_STRINGENCY "$VALIDATION_STRINGENCY" \
		--TMP_DIR "$TMPDIR"; then

		echo "ERROR: Duplicate marking failed for $bam_file"
		return
	fi

	echo "Successfully processed $bam_file"
}

active_jobs=0

for bam_file in "$OUTPUT_DIR"/*.bam; do
	[[ -e "$bam_file" ]] || continue

	mark_duplicates "$bam_file" &
	((active_jobs++))

	if [[ "$active_jobs" -ge "$MAX_JOBS" ]]; then
		wait -n
		((active_jobs--))
	fi
done

wait

echo "Duplicate marking completed for all BAM files in $OUTPUT_DIR"
