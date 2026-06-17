# Reference panel construction

This document describes how the 1000 Genomes Project reference panel used by the workflow was constructed. The benchmark analysis used a European-ancestry subset of the [1000 Genomes Project high-coverage GRCh38 phased ](https://www.internationalgenome.org/data-portal/data-collection/1000genomes_30x)  VCF files.

## 1. Download source 1000 Genomes Project files

Create a directory for the original 1000 Genomes Project VCF files:

    mkdir -p Data/reference_panel/1KG_30x_hg38/VCF
    cd Data/reference_panel/1KG_30x_hg38/VCF

Download chromosome-specific phased VCF files and their tabix indexes:

    BASE_URL="http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20220422_3202_phased_SNV_INDEL_SV"

    for chr in {1..22}; do
        wget -c "${BASE_URL}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
        wget -c "${BASE_URL}/1kGP_high_coverage_Illumina.chr${chr}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi"
    done

    wget -c "${BASE_URL}/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz"
    wget -c "${BASE_URL}/1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi"
    wget -c "${BASE_URL}/README_1kGP_phased_panel_110722.pdf"

The downloaded files should have the following structure:

    Data/reference_panel/1KG_30x_hg38/VCF/
    ├── 1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
    ├── 1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
    ├── ...
    ├── 1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
    ├── 1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel.vcf.gz.tbi
    ├── 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz
    └── 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2.vcf.gz.tbi

## 2. Sample metadata

The sample metadata file used to select ancestry-specific samples is provided in the repository:

    Data/example/1kg_samples_metadata.tsv

This file contains 1000 Genomes Project sample identifiers and population annotations. The reference panel construction script expects the following columns:

    column 1: sample ID
    column 2: sex
    column 6: superpopulation code


## 3. Build the European reference panel
For the benchmark analysis, samples with `EUR` in column 6 were selected. 

### <u>Users can adjust this part by changing the superpopulation code to one of the codes specified at 1000 Genomes documentation</u>

Run:

    scripts/build_1kgp_reference_panel.sh \
        --metadata Data/reference_panel/1KG_30x_hg38/1kg_samples_metadata.tsv \
        --vcf-dir Data/reference_panel/1KG_30x_hg38/VCF \
        --out-dir Data/reference_panel/EUR_nochr \
        --bref-jar /path/to/bref3.27Feb25.75f.jar \
        --superpopulation EUR

`bref3.27Feb25.75f.jar` is a part of Beagle 5.5 software which can downloaded from [Beagle website](https://faculty.washington.edu/browning/beagle/beagle.html).