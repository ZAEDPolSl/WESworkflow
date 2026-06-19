#!/bin/bash
set -euo pipefail

BAM_DIR="Data/example/original_bams"
FASTQ_DIR="Data/example/fastq"
TMP_DIR="Data/example/tmp"
SOURCE_BED="Data/bed/refGene_exons_splice5.nochr.bed"
REGIONS_PER_CHR=20

mkdir -p "$FASTQ_DIR" "$TMP_DIR"

if [[ ! -f "$SOURCE_BED" ]]; then
	echo "ERROR: BED file not found: $SOURCE_BED" >&2
	exit 1
fi

samples=(
	"sample1 WES_IL_N_1"
	"sample2 WES_IL_N_2"
	"sample3 WES_FD_N_1"
	"sample4 WES_FD_N_2"
)

echo "Creating example FASTQ files"
echo "BAM directory: $BAM_DIR"
echo "FASTQ directory: $FASTQ_DIR"
echo "Source BED: $SOURCE_BED"
echo "Regions per chromosome: $REGIONS_PER_CHR"
echo

for item in "${samples[@]}"; do
	out="$(echo "$item" | cut -d' ' -f1)"
	src="$(echo "$item" | cut -d' ' -f2)"
	bam="$BAM_DIR/${src}.bam"

	echo "Processing $out from $src"

	if [[ ! -f "$bam" ]]; then
		echo "ERROR: BAM file not found: $bam" >&2
		exit 1
	fi

	prefix="$(samtools idxstats "$bam" | awk 'BEGIN{p=""} $1=="chr1"{p="chr"} END{print p}')"
	extract_bed="$TMP_DIR/${out}.extract_regions.bed"

	awk -v n="$REGIONS_PER_CHR" '
		BEGIN {
			split("1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 X", chr_order)
			for (i in chr_order) wanted[chr_order[i]]=1
		}
		$1 in wanted && count[$1] < n {
			print $1, $2, $3
			count[$1]++
		}
	' OFS='\t' "$SOURCE_BED" > "$extract_bed"

	if [[ "$prefix" == "chr" ]]; then
		awk 'BEGIN{OFS="\t"} {$1="chr"$1; print}' "$extract_bed" > "$TMP_DIR/${out}.extract_regions.chr.bed"
		extract_bed="$TMP_DIR/${out}.extract_regions.chr.bed"
	fi

	echo "Selected regions: $(wc -l < "$extract_bed")"

	samtools view -F 2304 -L "$extract_bed" "$bam" | cut -f1 | sort -u > "$TMP_DIR/${out}.read_names.txt"
	echo "Selected read names: $(wc -l < "$TMP_DIR/${out}.read_names.txt")"

	if [[ ! -s "$TMP_DIR/${out}.read_names.txt" ]]; then
		echo "ERROR: No reads found for $out." >&2
		exit 1
	fi

	samtools view -b -F 2304 -N "$TMP_DIR/${out}.read_names.txt" "$bam" > "$TMP_DIR/${out}.subset.bam"
	echo "Subset BAM size: $(du -h "$TMP_DIR/${out}.subset.bam" | cut -f1)"

	samtools collate -@ 4 -u -O "$TMP_DIR/${out}.subset.bam" \
	| samtools fastq -@ 4 -1 "$FASTQ_DIR/${out}_1.fastq" -2 "$FASTQ_DIR/${out}_2.fastq" -0 /dev/null -s /dev/null -n -

	if [[ ! -s "$FASTQ_DIR/${out}_1.fastq" || ! -s "$FASTQ_DIR/${out}_2.fastq" ]]; then
		echo "ERROR: FASTQ files were not created correctly for $out." >&2
		exit 1
	fi

	gzip -f "$FASTQ_DIR/${out}_1.fastq" "$FASTQ_DIR/${out}_2.fastq"
	echo "Created:"
	ls -lh "$FASTQ_DIR/${out}_1.fastq.gz" "$FASTQ_DIR/${out}_2.fastq.gz"
	echo
done

echo "Done."
ls -lh "$FASTQ_DIR"
