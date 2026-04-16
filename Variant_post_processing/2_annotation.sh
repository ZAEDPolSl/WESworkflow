#!/bin/bash


# This script performs post-genotype-imputation VCF processing, with each chromosome being processed on another thread:
# (1) Annotates variants per chromosome using ANNOVAR and retains only functional variants (exonic/splicing)
# (2) Intersects imputed VCFs with functionally filtered original variants using bcftools isec
# (3) Annotates filtered VCFs with CADD v1.7 scores using bcftools annotate
#
# Input: per-chromosome VCF files (original + imputed)
# Output: functionally filtered and CADD-annotated VCF files, bgzip-compressed and indexed

WORKDIR="/path/to/results/dir" # all the downstream subdirectories with intermediate files will be created here
ANNOVAR_DB="/path/to/annovar/humandb"
mkdir -p $WORKDIR/Annotation
mkdir -p $WORKDIR/Intermediate


# step 3: Annotate variants in orignal VCF with ANNOVAR and filter them basing on function--------------------------------------------
INPUT_DIR="$WORKDIR/Imputation/Input_Chromosomes"
OUT_DIR="$WORKDIR/Annotation/multianno"
mkdir -p $OUT_DIR

export INPUT_DIR OUT_DIR ANNOVAR_DB
parallel -j 23 '
  table_annovar.pl \
    $INPUT_DIR/chr{1}.vcf.gz \
    $ANNOVAR_DB \
    -buildver hg38 \
    -out $OUT_DIR/chr{1} \
    -remove \
    --protocol refGene,exac03,gnomad41_exome,ALL.sites.2015_08,avsnp150,clinvar_20241215,dbnsfp47a \
    --operation g,f,f,f,f,f,f -nastring . -polish -vcfinput -thread 3

  grep "^#" $OUT_DIR/chr{1}.hg38_multianno.vcf > $OUT_DIR/chr{1}.hg38_multianno.function_filtered.vcf

  grep -v "^#" $OUT_DIR/chr{1}.hg38_multianno.vcf | \
  grep -E "Func.refGene=(exonic|splicing|exonic;splicing|ncRNA_exonic;splicing)" \
  >> $OUT_DIR/chr{1}.hg38_multianno.function_filtered.vcf

  bgzip -f $OUT_DIR/chr{1}.hg38_multianno.function_filtered.vcf
  tabix -f -p vcf $OUT_DIR/chr{1}.hg38_multianno.function_filtered.vcf.gz
' ::: {1..22} X



# step 4: Filter imputed file to the positions present in original vcf -----------------

export WORKDIR
parallel -j 23 '
  mkdir -p $WORKDIR/Intermediate/chr{1}_isec_tmp
  bcftools isec -Oz -p $WORKDIR/Intermediate/chr{1}_isec_tmp -n=2 -w1 \
    $WORKDIR/Imputation/Output_Chromosomes/Beagle_1kG_imputed/chr{1}.vcf.gz \
    $WORKDIR/Annotation/multianno/chr{1}.hg38_multianno.function_filtered.vcf.gz

  mv $WORKDIR/Intermediate/chr{1}_isec_tmp/0000.vcf.gz $WORKDIR/Imputation/Output_Chromosomes/Beagle_1kG_imputed/chr{1}.filtered.vcf.gz
  bcftools index $WORKDIR/Imputation/Output_Chromosomes/Beagle_1kG_imputed/chr{1}.filtered.vcf.gz
' ::: {1..22} X 




# step 5: add cadd scores to the orginial VCF ------------------------------------------------------
CADD_TSV="CADD_TSV="/path/to/CADD/whole_genome_SNVs.tsv.gz""
HEADER="$WORKDIR/Annotation/cadd.header.txt"
IN_DIR="$WORKDIR/Annotation/multianno"
OUT_DIR="$WORKDIR/Annotation/CADD_scored"

# header INFO for CADD
cat > $WORKDIR/Annotation/cadd.header.txt <<'EOF'
##INFO=<ID=CADD1.7_RAW,Number=1,Type=Float,Description="CADD RawScore">
##INFO=<ID=CADD1.7_PHRED,Number=1,Type=Float,Description="CADD PHRED">
EOF

mkdir -p "$OUT_DIR"

export IN_DIR OUT_DIR CADD_TSV HEADER
parallel -j 23 '
  IN_VCF=$IN_DIR/chr{1}.hg38_multianno.function_filtered.vcf.gz
  OUT_VCF=$OUT_DIR/chr{1}.hg38_multianno.filtered_CADD_scored.vcf.gz

  bcftools annotate \
    -a $CADD_TSV \
    --threads 3 \
    -c CHROM,POS,REF,ALT,INFO/CADD1.7_RAW:=5,INFO/CADD1.7_PHRED:=6 \
    -h $HEADER \
    -O z -o $OUT_VCF \
    $IN_VCF && \
  bcftools index -f $OUT_VCF
' ::: {1..22} X


