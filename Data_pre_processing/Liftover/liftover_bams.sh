#!/bin/bash

# This script performs liftover of BAM files from hg19 to hg38 using CrossMap tool and GNU parallel for batch processing
#/general/dataset/directory/    <- $WORKDIR
#├── fastq_files/           
#├── trimmed_fastq_files/      
#│   └── fastqc/ 
#├── bam_files/
#├── hg19_bam_files/  

export TMPDIR="optional/tmp/dir" # set a temporary directory for CrossMap
mkdir -p "$TMPDIR"

THREADS=24 # number of threads for parallel computing
CHAIN="/path/to/chain_file/GRCh37_to_GRCh38.chain.gz"
INPUT_DIR="/path/to/hg19_bam_files" # path to hg19 bams
OUTPUT_DIR="/path/for/lifted/bam_files" # path to output hg38 bams

mkdir -p "$OUTPUT_DIR" # create output directory

process_sample() {
    line="$1"
    sample=$(basename "$line" .bam) 

    out_bam="$OUTPUT_DIR/${sample}"

    echo "Processing $sample …"

    if [ -f "${out_bam}.bam" ]; then
        echo "   -> exists, skipping"
        return
    fi

    CrossMap bam --chromid s "$CHAIN" "$line" "$out_bam"
}

export -f process_sample
export CHAIN OUTPUT_DIR TMPDIR


find "$INPUT_DIR" -type f -depth 2 -name "*.bam" | parallel --will-cite -j $THREADS 'bash -lc "process_sample {}"'