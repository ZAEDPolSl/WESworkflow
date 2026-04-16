#!/bin/bash
# This script is used to run GLnexus for merging and genotyping gVCF files

WORKDIR="path/to/results/dir" # all the downstream subdirectories with intermediate files will be created here
OUTPUT_DIR="$WORKDIR/Genotyping"
REFERENCE="/path/to/the/hg38_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"

# path to the regions of restriction
REGIONS="Data/bed/refGene_exons_splice5.nochr.bed"

SCRATCH_DIR="optional/path/to/tmp/GLnexus_scratch" # path to GLnexus intermediate files

THREADS=40

mkdir -p $OUTPUT_DIR

rm -rf $SCRATCH_DIR

# catch all the snv files in the project
readarray -t INPUT_DIRS < <(find /path/to/general/datasets/parent/directory -maxdepth 5 -type d -name "snv")


FILES=$(find "${INPUT_DIRS[@]}" -type f -name "*.g.vcf.gz")

# Joint genotyping
glnexus_cli \
    --config DeepVariantWES \
    --bed "$REGIONS" \
    --dir "$SCRATCH_DIR" \
    $FILES \
    > "$OUTPUT_DIR/genotyped_multisample.bcf"
 
# bcf to vcf conversion
bcftools view \
    --threads $THREADS \
    -Ov $OUTPUT_DIR/genotyped_multisample.bcf \
    -o $OUTPUT_DIR/genotyped_multisample.vcf

bcftools view \
--threads $THREADS \
-Oz $OUTPUT_DIR/genotyped_multisample.bcf \
-o $OUTPUT_DIR/genotyped_multisample.vcf.gz

# index the gzipped vcf
bcftools index --threads 3 -f $OUTPUT_DIR/genotyped_multisample.vcf.gz