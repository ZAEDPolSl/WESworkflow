#!/usr/bin/env bash
set -euo pipefail

# This script performs liftover of BAM files from GRCh37/hg19 to GRCh38
# using CrossMap and GNU parallel for batch processing.
# Usage: bash Data_pre_processing/Liftover/liftover_bams.sh config/local_config.yaml

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

TMPDIR="$(python "$READ_CONFIG" "$CONFIG" directories.temporary_dir)"
THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.liftover_threads)"
CHAIN="$(python "$READ_CONFIG" "$CONFIG" resources.liftover_chain)"
INPUT_DIR="$(python "$READ_CONFIG" "$CONFIG" directories.hg19_bam)"
OUTPUT_DIR="$(python "$READ_CONFIG" "$CONFIG" directories.lifted_bam)"

CHROMID="s"

export TMPDIR

if [[ ! -d "$INPUT_DIR" ]]; then
	echo "ERROR: Input BAM directory not found: $INPUT_DIR"
	exit 1
fi

if [[ ! -f "$CHAIN" ]]; then
	echo "ERROR: Liftover chain file not found: $CHAIN"
	exit 1
fi

mkdir -p "$OUTPUT_DIR" "$TMPDIR"

process_sample() {
	local bam_file="$1"
	local sample
	local out_bam

	sample="$(basename "$bam_file" .bam)"
	out_bam="$OUTPUT_DIR/${sample}.bam"

	echo "Processing $sample"

	if [[ -f "$out_bam" ]]; then
		echo "Skipping $sample - output exists"
		return
	fi

	CrossMap bam --chromid "$CHROMID" "$CHAIN" "$bam_file" "$out_bam"
}

export -f process_sample
export CHAIN OUTPUT_DIR TMPDIR CHROMID

bam_count="$(find "$INPUT_DIR" -type f -name "*.bam" | wc -l)"

if [[ "$bam_count" -eq 0 ]]; then
	echo "ERROR: No BAM files found in $INPUT_DIR"
	exit 1
fi

find "$INPUT_DIR" -type f -name "*.bam" -print0 | \
	parallel --will-cite -0 -j "$THREADS" bash -c 'process_sample "$1"' _ {}
