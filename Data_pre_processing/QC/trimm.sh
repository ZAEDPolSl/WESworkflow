#!/usr/bin/env bash
set -euo pipefail

# Trims paired-end FASTQ files with Trimmomatic, then runs FastQC and MultiQC.
# Usage: bash Data_pre_processing/QC/trimm.sh config/local_config.yaml
# Required config fields: directories.raw_fastq, directories.trimmed_fastq,
# resources.trimmomatic_adapters, parameters.trimmomatic_* and parameters.fastqc_threads.

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

INPUT_DIR="$(python "$READ_CONFIG" "$CONFIG" directories.raw_fastq)"
OUTPUT_DIR="$(python "$READ_CONFIG" "$CONFIG" directories.trimmed_fastq)"
FASTQC_DIR="$OUTPUT_DIR/fastqc"

TRIMMO_THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_threads)"
FASTQC_THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.fastqc_threads)"
ADAPTERS="$(python "$READ_CONFIG" "$CONFIG" resources.trimmomatic_adapters)"
ILLUMINACLIP="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_illuminaclip)"
SLIDINGWINDOW="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_slidingwindow)"
MINLEN="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_minlen)"

if [[ ! -d "$INPUT_DIR" ]]; then
	echo "ERROR: Input FASTQ directory not found: $INPUT_DIR"
	exit 1
fi

if [[ ! -f "$ADAPTERS" ]]; then
	echo "ERROR: Trimmomatic adapter file not found: $ADAPTERS"
	exit 1
fi

mkdir -p "$OUTPUT_DIR" "$FASTQC_DIR"

i=1
total="$(find "$INPUT_DIR" -name "*_1.fastq.gz" | wc -l)"

if [[ "$total" -eq 0 ]]; then
	echo "ERROR: No *_1.fastq.gz files found in $INPUT_DIR"
	exit 1
fi

find "$INPUT_DIR" -name "*_1.fastq.gz" | sort | while read -r r1; do
	sample="$(basename "$r1" _1.fastq.gz)"
	r2="$INPUT_DIR/${sample}_2.fastq.gz"

	if [[ ! -f "$r2" ]]; then
		echo "WARNING: Missing R2 file for sample $sample: $r2"
		continue
	fi

	echo "Processing sample $sample ($i/$total)"

	output_file="$OUTPUT_DIR/${sample}_1.fastq.gz"
	if [[ -f "$output_file" ]]; then
		echo "Skipping $sample - output exists"
		((i++))
		continue
	fi

	trimmomatic PE -threads "$TRIMMO_THREADS" \
		"$r1" "$r2" \
		"$OUTPUT_DIR/${sample}_1.fastq.gz" "$OUTPUT_DIR/${sample}_1_unpaired.fastq.gz" \
		"$OUTPUT_DIR/${sample}_2.fastq.gz" "$OUTPUT_DIR/${sample}_2_unpaired.fastq.gz" \
		"ILLUMINACLIP:${ADAPTERS}:${ILLUMINACLIP}" \
		"SLIDINGWINDOW:${SLIDINGWINDOW}" "MINLEN:${MINLEN}"

	((i++))
done

fastqc --threads "$FASTQC_THREADS" --outdir "$FASTQC_DIR" "$OUTPUT_DIR"/*[12].fastq.gz

cd "$FASTQC_DIR"
multiqc .