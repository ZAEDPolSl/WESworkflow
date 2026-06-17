#!/usr/bin/env bash
set -euo pipefail

META=""
VCF_DIR=""
OUT_DIR=""
BREF_JAR=""
SUPERPOP="EUR"
JOBS=22
THREADS=10
XMX="12g"
MAC=5

usage() {
	cat <<'EOF'
Usage:
  scripts/build_1kgp_reference_panel.sh \
    --metadata Data/example/1kg_samples_metadata.tsv \
    --vcf-dir Data/reference_panel/1KG_30x_hg38/VCF \
    --out-dir Data/reference_panel/EUR_nochr \
    --bref-jar /path/to/bref3.jar \
    [--superpopulation EUR] [--jobs 22] [--threads 10] [--xmx 12g] [--mac 5]

Metadata format:
  column 1: sample ID
  column 2: sex
  column 6: superpopulation code
EOF
}

while [[ $# -gt 0 ]]; do
	case "$1" in
		--metadata) META="$2"; shift 2 ;;
		--vcf-dir) VCF_DIR="$2"; shift 2 ;;
		--out-dir) OUT_DIR="$2"; shift 2 ;;
		--bref-jar) BREF_JAR="$2"; shift 2 ;;
		--superpopulation) SUPERPOP="$2"; shift 2 ;;
		--jobs) JOBS="$2"; shift 2 ;;
		--threads) THREADS="$2"; shift 2 ;;
		--xmx) XMX="$2"; shift 2 ;;
		--mac) MAC="$2"; shift 2 ;;
		-h|--help) usage; exit 0 ;;
		*) echo "Unknown argument: $1"; usage; exit 1 ;;
	esac
done

if [[ -z "$META" || -z "$VCF_DIR" || -z "$OUT_DIR" || -z "$BREF_JAR" ]]; then
	echo "Missing required argument."
	usage
	exit 1
fi

[[ -f "$META" ]] || { echo "Metadata file not found: $META"; exit 1; }
[[ -d "$VCF_DIR" ]] || { echo "VCF directory not found: $VCF_DIR"; exit 1; }
[[ -f "$BREF_JAR" ]] || { echo "BREF3 JAR not found: $BREF_JAR"; exit 1; }

BREF_DIR="$OUT_DIR/BREF3"
mkdir -p "$OUT_DIR" "$BREF_DIR"

SAMPLES="$OUT_DIR/${SUPERPOP}.samples.txt"
FEMALE_SAMPLES="$OUT_DIR/${SUPERPOP}.female.samples.txt"
CHRMAP="$OUT_DIR/chr_rename.txt"

awk -v sp="$SUPERPOP" -F'\t' '$6==sp{print $1}' "$META" > "$SAMPLES"
awk -v sp="$SUPERPOP" -F'\t' '$6==sp && $2=="female"{print $1}' "$META" > "$FEMALE_SAMPLES"

[[ -s "$SAMPLES" ]] || { echo "No samples found for superpopulation: $SUPERPOP"; exit 1; }
[[ -s "$FEMALE_SAMPLES" ]] || echo "Warning: no female samples found for chromosome X."

{ for i in {1..22}; do echo -e "chr${i}\t${i}"; done; echo -e "chrX\tX"; } > "$CHRMAP"

export VCF_DIR OUT_DIR BREF_DIR BREF_JAR SAMPLES CHRMAP SUPERPOP THREADS XMX MAC

parallel -j "$JOBS" --linebuffer '
	chr={1}
	in_vcf="$VCF_DIR/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
	out_vcf="$OUT_DIR/1kGP.${chr}.${SUPERPOP}.nochr.vcf.gz"
	out_bref="$BREF_DIR/1kGP.${chr}.${SUPERPOP}.nochr.bref3"

	[[ -f "$in_vcf" ]] || { echo "Input VCF not found: $in_vcf"; exit 1; }

	bcftools view --threads "$THREADS" -M2 -m2 -v snps,indels -i "MAC>=${MAC}" -S "$SAMPLES" "$in_vcf" -Ou \
		| bcftools annotate --rename-chrs "$CHRMAP" -Oz -o "$out_vcf"

	tabix -f -p vcf "$out_vcf"
	zcat "$out_vcf" | java -Xmx"$XMX" -jar "$BREF_JAR" > "$out_bref"
' ::: {1..22}

X_VCF="$VCF_DIR/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"

if [[ -f "$X_VCF" && -s "$FEMALE_SAMPLES" ]]; then
	out_vcf="$OUT_DIR/1kGP.X.${SUPERPOP}.nochr.vcf.gz"
	out_bref="$BREF_DIR/1kGP.X.${SUPERPOP}.nochr.bref3"

	bcftools view --threads "$THREADS" -S "$FEMALE_SAMPLES" --force-samples -r chrX "$X_VCF" -Ou \
		| bcftools annotate --rename-chrs "$CHRMAP" -Oz -o "$out_vcf"

	tabix -f -p vcf "$out_vcf"
	zcat "$out_vcf" | java -Xmx"$XMX" -jar "$BREF_JAR" > "$out_bref"
else
	echo "Skipping chromosome X: input VCF or female sample list not available."
fi

echo "Reference panel written to: $OUT_DIR"
