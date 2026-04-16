#!/bin/bash

# This script performs trimming of paired-end FASTQ files using Trimmomatic. 
# It is to be used only if initial QC does not pass.
# It also runs FastQC on the trimmed files and generates a MultiQC report in the end.

#/general/dataset/directory/    <- $WORKDIR
#├── fastq_files/           
#├── trimmed_fastq_files/      
#│   └── fastqc/ 

TRIMMO_THREADS=8
FASTQC_THREADS=24
WORKDIR="path/to/dataset/dir" 
INPUT_DIR=$WORKDIR/fastq_files
OUTPUT_DIR=$WORKDIR/trimmed_fastq_files
FASTQC_DIR=$OUTPUT_DIR/fastqc

mkdir -p $OUTPUT_DIR $FASTQC_DIR

i=1
total=$(find "$INPUT_DIR" -name "*_1.fastq.gz" | wc -l) # total samples to process

# ADJUST FILE EXTENSION
find "$INPUT_DIR" -name "*_1.fastq.gz" | while read -r r1; do

    sample=$(basename "$r1" _1.fastq.gz)
    r2="$INPUT_DIR/${sample}_2.fastq.gz"


    echo "Processing sample $sample ($i/$total)"
    # skip if sample processed
    output_file="$OUTPUT_DIR/${sample}_1.fastq"
    if [ -f "$output_file" ]; then
        echo "Skipping $sample – output exists"
        ((i++))
        continue
    fi

# quality and adapters trimming, adjust if needed
    trimmomatic PE -threads $TRIMMO_THREADS \
        "$r1" "$r2" \
        "$OUTPUT_DIR/${sample}_1.fastq.gz" "$OUTPUT_DIR/${sample}_1_unpaired.fastq.gz" \
        "$OUTPUT_DIR/${sample}_2.fastq.gz" "$OUTPUT_DIR/${sample}_2_unpaired.fastq.gz" \
        ILLUMINACLIP:$CONDA_PREFIX/share/trimmomatic-0.40-0/adapters/TruSeq3-PE.fa:2:30:10 \ 
        SLIDINGWINDOW:4:20 MINLEN:45 
        
    ((i++))
done

# run quality check
fastqc \
    --threads $FASTQC_THREADS \
    --outdir "$FASTQC_DIR" \
    "$OUTPUT_DIR"/*[12].fastq.gz

# move to the FastQC output directory and run MultiQC to generate a report
cd $FASTQC_DIR
multiqc .
