#!/bin/bash
set -euo pipefail

OUT_DIR="Data/example/original_bams"
BASE_URL="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/seqc/Somatic_Mutation_WG/data/WES"

mkdir -p "$OUT_DIR"

samples=(
	"WES_IL_N_1"
	"WES_IL_N_2"
	"WES_FD_N_1"
	"WES_FD_N_2"
)

for sample in "${samples[@]}"; do
	wget -c -O "$OUT_DIR/${sample}.bam" "$BASE_URL/${sample}.bwa.dedup.bam"
	wget -c -O "$OUT_DIR/${sample}.bam.bai" "$BASE_URL/${sample}.bwa.dedup.bai"
done
