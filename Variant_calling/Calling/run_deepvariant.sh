#!/bin/bash

# This script performs variant calling on Whole Exome Sequencing BAMs using docker image DeepVariant tool 
# producing VCF and GVCF files
#/general/dataset/directory/    <- $WORKDIR
#├── fastq_files/           
#├── trimmed_fastq_files/      
#│   └── fastqc/ 
#├── bam_files/    
#├── vcf_files/      
#│   └── snv/ 

WORKDIR="path/to/dataset/dir"
INPUT_DIRS="$WORKDIR/bam_files"
OUTPUT_DIR="$WORKDIR/vcf_files/snv"

REFERENCE_DIR="/path/to/the/hg38_reference"
REFERENCE_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# path to the regions of restriction
REGIONS_DIR="Data/bed"
REGIONS_FILE="refGene_exons_splice5.nochr.bed"

mkdir -p "$OUTPUT_DIR" "$OUTPUT_DIR/logs"

ERROR_LOG="$OUTPUT_DIR/logs/error.log"

THREADS=30

BIN_VERSION="1.8.0"


# CHECK HOW YOUR BAMS ARE NAMED BECAUSE SCRIPT MATCHES THE SUFFIX
total=$(find "$INPUT_DIRS" -name "*sorted.bam" | wc -l)
i=0

find $INPUT_DIRS -name "*sorted.bam" |
while IFS= read -r line; do
    filename=$(basename "$line")
    sample_name="${filename%sorted.bam}"
    ((i++))
    echo "[$i/$total] Processing $sample_name …"

    if [ -f "$OUTPUT_DIR/${sample_name}.vcf.gz" ]; then
        echo "    → exists, skipping"
        continue
    fi

    docker run --rm --gpus 1 \
    -v "${INPUT_DIRS}":"/input" \
    -v "${OUTPUT_DIR}":"/output" \
    -v "${REFERENCE_DIR}":"/reference" \
    -v "${REGIONS_DIR}":"/regions" \
    -v /mnt/hot/tmp:/tmp \
    google/deepvariant:"${BIN_VERSION}-gpu" \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type WES \
    --ref "/reference/$REFERENCE_FILE" \
    --reads "/input/$filename" \
    --regions "/regions/$REGIONS_FILE" \
    --output_vcf "/output/${sample_name}.vcf.gz" \
    --output_gvcf "/output/${sample_name}.g.vcf.gz" \
    --logging_dir="/output/logs" \
    --num_shards "$THREADS" \
    --intermediate_results_dir "/tmp/dv_${sample_name}" \
    --sample_name="$sample_name" 
done 






