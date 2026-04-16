#!/bin/bash

# This script performs alignment of paired-end FASTQ files to a reference genome using BWA MEM and sorts the output BAM files.
# When the alignment is finished, it performs duplicates marking

#/general/dataset/directory/    <- $WORKDIR
#├── fastq_files/           
#├── trimmed_fastq_files/      
#│   └── fastqc/ 
#├── bam_files/    

THREADS=60 # number of threads to use for alignment and sorting

REFERENCE_DIR="/path/to/the/hg38_reference"
REFERENCE_FILE="Homo_sapiens.GRCh38.dna.primary_assembly.fa"

WORKDIR="path/to/dataset/dir" 
INPUT_DIR=$WORKDIR/fastq_files
OUTPUT_DIR=$WORKDIR/bam_files


TMPDIR="optional/tmp/dir" # Temporary directory for sorting
mkdir -p "$TMPDIR"
export TMPDIR


mkdir -p $OUTPUT_DIR


echo "STARTED PROCESSING"

# store failed samples
failed_samples=()
i=1
total=$(find "$INPUT_DIR" -name "*_1.fastq" | wc -l) # total samples to process



for r1 in $INPUT_DIR/*_1.fastq; do
    sample=$(basename "$r1" _1.fastq) 
    # sample=${sample:0:16} uncomment if TCGA, take only 1st 16 letters!
    r2="$INPUT_DIR/${sample}_2.fastq"

    # Check if the BAM file already exists skip if yes
    if [ -f "$OUTPUT_DIR/${sample}.bam" ]; then
        ((i++))
        continue
    fi

    echo "Processing sample $sample ($i/$total)"

    # Align and sort the BAM file 
    bwa mem -t $THREADS -T 0 -v 1 \
        -R "@RG\tID:${sample}\tSM:${sample}\tLB:WES\tPL:Illumina" \
        "$REFERENCE_DIR/$REFERENCE_FILE" "$r1" "$r2" | \
    samtools view -Shb -@ $THREADS | \
    samtools sort -@ $THREADS -T "$TMPDIR/${sample}_tmp" -o "$OUTPUT_DIR/${sample}.bam"

    

    if [ $? -ne 0 ]; then
        echo "Error during alignment or sorting for sample $sample"
        failed_samples+=("$sample")
        ((i++))
        continue
    fi


    ((i++))
done

# Summary of failed samples
if [ ${#failed_samples[@]} -gt 0 ]; then
    echo "The following samples failed to process:"
    for sample in "${failed_samples[@]}"; do
        echo "- $sample"
    done
else
    echo "ALIGNMENT COMPLETED SUCCESSFULLY"
fi

# MARK DUPLICATES ======================================
MAX_JOBS=50  

mark_duplicates() {
    local bam_file="$1"
    local output_bam="${bam_file%.bam}_marked.bam"
    local metrics_file="${bam_file%.bam}_metrics.txt"

    echo "Processing $bam_file"
    picard MarkDuplicates \
        --INPUT "$bam_file" \
        --OUTPUT "$output_bam" \
        --METRICS_FILE "$metrics_file" \
        --CREATE_INDEX true \
        --VALIDATION_STRINGENCY STRICT \
        --TMP_DIR "$TMPDIR"

    if [ $? -ne 0 ]; then
        echo "Error marking duplicates for $bam_file"
    else
        echo "Successfully processed $bam_file"
    fi

    rm "$metrics_file" 
}


active_jobs=0


for bam_file in "$OUTPUT_DIR"/*.bam; do

    output_bam="${bam_file%.bam}_marked.bam"

    if [[ -e "$output_bam" || "$bam_file" == *_marked.bam ]]; then
        echo "Skipping $bam_file"
        continue
    fi

    mark_duplicates "$bam_file" &  
    ((active_jobs++))

    if [ "$active_jobs" -ge "$MAX_JOBS" ]; then
        wait -n   
        ((active_jobs--))
    fi
done

echo "Marking duplicates completed for all BAM files in $OUTPUT_DIR."

rm -rf "$TMPDIR"
