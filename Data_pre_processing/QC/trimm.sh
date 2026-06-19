#!/usr/bin/env bash
set -euo pipefail

CONFIG="${1:-config/local_config.yaml}"
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"

if [[ "$CONFIG" != /* ]]; then
	CONFIG="$REPO_DIR/$CONFIG"
fi

READ_CONFIG="$REPO_DIR/scripts/read_config.py"

if [[ ! -f "$CONFIG" ]]; then
	echo "ERROR: Config file not found: $CONFIG" >&2
	exit 1
fi

resolve_path() {
	local path="$1"

	if [[ "$path" == /* ]]; then
		echo "$path"
	else
		echo "$REPO_DIR/$path"
	fi
}

INPUT_DIR="$(resolve_path "$(python "$READ_CONFIG" "$CONFIG" directories.raw_fastq)")"
OUTPUT_DIR="$(resolve_path "$(python "$READ_CONFIG" "$CONFIG" directories.trimmed_fastq)")"
FASTQC_DIR="$OUTPUT_DIR/fastqc"

TRIMMO_THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_threads)"
FASTQC_THREADS="$(python "$READ_CONFIG" "$CONFIG" parameters.fastqc_threads)"
TRIM_ADAPTERS="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_trim_adapters 2>/dev/null || echo true)"
ADAPTER_FILE="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_adapter_file 2>/dev/null || echo TruSeq3-PE.fa)"
ILLUMINACLIP="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_illuminaclip)"
SLIDINGWINDOW="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_slidingwindow)"
MINLEN="$(python "$READ_CONFIG" "$CONFIG" parameters.trimmomatic_minlen)"

TRIM_ADAPTERS="$(echo "$TRIM_ADAPTERS" | tr '[:upper:]' '[:lower:]')"

if [[ ! -d "$INPUT_DIR" ]]; then
	echo "ERROR: Input FASTQ directory not found: $INPUT_DIR" >&2
	exit 1
fi

mkdir -p "$OUTPUT_DIR" "$FASTQC_DIR"

TRIMMOMATIC_STEPS=()

if [[ "$TRIM_ADAPTERS" == "true" ]]; then
	if [[ -z "${CONDA_PREFIX:-}" ]]; then
		echo "ERROR: CONDA_PREFIX is not set. Activate the Conda environment first." >&2
		exit 1
	fi

	ADAPTER_PATH="$CONDA_PREFIX/share/trimmomatic/adapters/$ADAPTER_FILE"

	if [[ ! -f "$ADAPTER_PATH" ]]; then
		echo "ERROR: Adapter file not found: $ADAPTER_PATH" >&2
		echo "Check Trimmomatic installation or parameters.trimmomatic_adapter_file." >&2
		exit 1
	fi

	TRIMMOMATIC_STEPS+=("ILLUMINACLIP:${ADAPTER_PATH}:${ILLUMINACLIP}")
elif [[ "$TRIM_ADAPTERS" != "false" ]]; then
	echo "ERROR: parameters.trimmomatic_trim_adapters must be true or false" >&2
	exit 1
fi

TRIMMOMATIC_STEPS+=("SLIDINGWINDOW:${SLIDINGWINDOW}" "MINLEN:${MINLEN}")

i=1
total="$(find "$INPUT_DIR" -name "*_1.fastq.gz" | wc -l)"

if [[ "$total" -eq 0 ]]; then
	echo "ERROR: No *_1.fastq.gz files found in $INPUT_DIR" >&2
	exit 1
fi

find "$INPUT_DIR" -name "*_1.fastq.gz" | sort | while read -r r1; do
	sample="$(basename "$r1" _1.fastq.gz)"
	r2="$INPUT_DIR/${sample}_2.fastq.gz"

	if [[ ! -f "$r2" ]]; then
		echo "WARNING: Missing R2 file for sample $sample: $r2" >&2
		continue
	fi

	echo "Processing sample $sample ($i/$total)"

	if [[ -f "$OUTPUT_DIR/${sample}_1.fastq.gz" && -f "$OUTPUT_DIR/${sample}_2.fastq.gz" ]]; then
		echo "Skipping $sample - output exists"
		((i++))
		continue
	fi

	trimmomatic PE -threads "$TRIMMO_THREADS" \
		"$r1" "$r2" \
		"$OUTPUT_DIR/${sample}_1.fastq.gz" "$OUTPUT_DIR/${sample}_1_unpaired.fastq.gz" \
		"$OUTPUT_DIR/${sample}_2.fastq.gz" "$OUTPUT_DIR/${sample}_2_unpaired.fastq.gz" \
		"${TRIMMOMATIC_STEPS[@]}"

	((i++))
done

fastqc --threads "$FASTQC_THREADS" --outdir "$FASTQC_DIR" "$OUTPUT_DIR"/*[12].fastq.gz

cd "$FASTQC_DIR"
multiqc .