#!/bin/bash

# This script performs genotype preprocessing and imputation:
# (1) Normalizes and filters a multisample VCF (split multiallelics, left-align, keep SNVs with AD > 5),
#     and splits data into per-chromosome VCF files
# (2) Conforms genotypes to a reference panel using conform-gt
# (3) Performs genotype imputation with Beagle using a population reference panel and genetic maps
#
# Input: multisample VCF (e.g. GLnexus output)
# Output: per-chromosome imputed and indexed VCF files
#
# Reference panel is available at:
# https://doi.org/10.5281/zenodo.19626919 and should be downloaded into "Data/reference_panel/EUR_nochr" path

# Genetic maps are not included and should be downloaded separately.

WORKDIR="path/to/results/dir" # all the downstream subdirectories with intermediate files will be created here
IN_VCF="$WORKDIR/Genotyping/genotyped_multisample.vcf.gz" # GLnexus output
REFERENCE="/path/to/the/hg38_reference/Homo_sapiens.GRCh38.dna.primary_assembly.fa"


mkdir -p $WORKDIR/Imputation
mkdir -p $WORKDIR/Imputation/Input_Chromosomes
mkdir -p $WORKDIR/Imputation/Output_Chromosomes


# step 1: Split multiallelic sites, left-align filter the records where any AD > 5 ans SNVs only, split into chromosomes --------------------------------------------

export WORKDIR REFERENCE IN_VCF
parallel -j 23 '
  bcftools norm -m -both -f ${REFERENCE} -d all -r {1} --threads 3 -Oz "$IN_VCF" | \
  bcftools view -v snps -i "MAX(AD[*]) > 5" --threads 3 -Oz \
    -o "$WORKDIR/Imputation/Input_Chromosomes/chr{1}.vcf.gz" && \
  bcftools index --threads 3 -f "$WORKDIR/Imputation/Input_Chromosomes/chr{1}.vcf.gz"
' ::: {1..22} X

echo "Preprocessing completed!"

# step 2: conform genotypes and impute with beagle ----------------------------------------------------------------------
BEAGLE_JAR="/path/to/beagle/beagle.jar"
CONFORM_JAR="/path/to/beagle/conform-gt.jar"
MAP_DIR="/path/to/beagle/genetic_maps"
INPUT_DIR="$WORKDIR/Imputation/Input_Chromosomes"
OUT_DIR="$WORKDIR/Imputation/Output_Chromosomes/Beagle_1kG_imputed"
REF_DIR="Data/reference_panel/EUR_nochr"

mkdir -p "$OUT_DIR"

 
export BEAGLE_JAR CONFORM_JAR MAP_DIR INPUT_DIR OUT_DIR REF_DIR

echo "Started genotype conformation..."

parallel --line-buffer -j 23 "
    echo START chr{1}
    
    java -jar \$CONFORM_JAR \
        gt=\$INPUT_DIR/chr{1}.vcf.gz \
        ref=\$REF_DIR/1kGP.{1}.EUR.nochr.vcf.gz \
        chrom={1} \
        match=POS \
    out=\$INPUT_DIR/chr{1}.conformed
    
    bcftools index -f "$INPUT_DIR/chr{1}.conformed.vcf.gz"

    echo DONE chr{1}
" ::: {1..22} X

echo "Done!"

parallel -j 23 "
  java -Xmx50g -jar \$BEAGLE_JAR \
    gt=\$INPUT_DIR/chr{1}.conformed.vcf.gz \
    map=\$MAP_DIR/plink.chr{1}.GRCh38.map \
    ref=\$REF_DIR/BREF3/1kGP.{1}.EUR.nochr.bref3 \
    nthreads=10 \
    out=\$OUT_DIR/chr{1}

    bcftools index -f "$OUT_DIR/chr{1}.vcf.gz"
" ::: {1..22} X