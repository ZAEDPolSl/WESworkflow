#!/bin/bash

# This script was used for generating the BED file with regions containig protein coding genes used in the project

# wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/refGene.txt.gz
# gunzip refGene.txt.gz

# convert to BED format
 awk -F'\t' 'BEGIN{OFS="\t"} {
    split($10, starts, ",")
    split($11, ends, ",")
    for (i = 1; i <= $9; i++) {
        print $3, starts[i], ends[i]
    }
}' refGene.txt | sort -k1,1 -k2,2n | bedtools merge > refGene_exons.bed

# clean up alternate contigs
egrep '^chr([1-9]|1[0-9]|2[0-2]|X|Y|M)[[:space:]]' refGene_exons.bed > tmp && mv tmp refGene_exons.bed

# add splicings:   +- 5 bases around exon/intron border
bedtools slop -b 5 -i refGene_exons.bed -g genome_file.txt | bedtools merge > refGene_exons_splice5.bed

# remove 'chr' prefix 
sed 's/^chr//' refGene_exons_splice5.bed > refGene_exons_splice5.nochr.bed

# show final chromosomes
cut -f1 refGene_exons_splice5.nochr.bed | sort | uniq

# change MT to M for GLnexus (may not be useful anymore)
awk 'BEGIN{OFS="\t"} $1 ~ /^(M|[1-9][0-9]*|X|Y)$/ { if ($1 == "M") $1 = "MT"; print }' refGene_exons_splice5.nochr.bed > GLnexus.bed
