# A practical workflow for correcting kit-specific effects in whole-exome sequencing data

**WESworkflow** is a modular pipeline for processing whole-exome sequencing (WES) data from raw reads or pre-aligned BAM files to gene-level variant feature matrices. The workflow performs quality control, trimming, alignment or liftover, duplicate marking, variant calling, joint genotyping, genotype imputation, functional annotation, CADD-weighted allele fraction (CWAF) aggregation, batch-effect assessment, and MNAR-aware gene-level imputation.

The workflow was developed for GRCh38/hg38 WES data and is intended for multi-source cohorts where exome capture kit differences may introduce systematic gene-level detection biases.

## Quick start

A detailed setup guide is provided in [**installation guide**](docs/installation.md).


In a typical setup, the user should:

```bash
git clone https://github.com/ZAEDPolSl/WESworkflow.git
cd WESworkflow
cp config/example_config.yaml config/local_config.yaml
```

Then edit `config/local_config.yaml` so that all tool paths, resource paths, working directories, and adjustable parameters match the local environment.

Before running the full workflow on a cohort, run the lightweight installation check or dry-run procedure described in [**installation guide**](docs/installation.md). This step is intended to verify that dependencies, external binaries, reference files, and configuration paths are visible before launching computationally expensive jobs.

## Configuration

The workflow uses a YAML configuration file:

```text
config/example_config.yaml
```

Copy this file to a local, user-specific configuration file:

```text
config/local_config.yaml
```

The configuration file contains:

- paths to external tools, such as Beagle, conform-gt, ANNOVAR, and the DeepVariant Docker image;
- paths to external resources, such as the reference genome, liftover chain file, Beagle genetic maps, reference panel, RefSeq BED file, and CADD prescored database;
- user-adjustable workflow parameters, including Trimmomatic settings, BWA alignment sensitivity, thread counts, job counts, DeepVariant runtime settings, Beagle settings, and input/output directories.

In the configuration file, `jobs` denotes externally parallelized processes, while `threads` denotes CPU threads used internally by a given tool.

Several method-specific analytical constants used after gene-level feature generation are kept in the corresponding R scripts rather than in the global YAML file. These are placed near the beginning of the relevant scripts to keep them visible while avoiding overloading the main configuration file.

## Requirements

The workflow requires a Unix-like environment with Bash, Conda, Docker or Singularity-compatible DeepVariant execution, Python, R, and several bioinformatics tools. The exact installation procedure and package installation commands are described in [**installation guide**](docs/installation.md).

Core command-line tools:

```text
FastQC
MultiQC
Trimmomatic
BWA
Samtools
Picard
CrossMap, only for GRCh37/hg19 to GRCh38/hg38 BAM liftover
DeepVariant
GLnexus
bcftools
htslib/bgzip/tabix
Beagle
conform-gt
ANNOVAR
GNU parallel
```

Python is used for configuration parsing, PARC clustering, and variant-to-gene feature calculation. R is used for feature loading, UMAP visualization, batch-effect metrics, Gaussian mixture modeling, and gene-level feature imputation.

Important Python/R packages used in downstream analyses include:

```text
Python: pandas, numpy, scikit-learn, cyvcf2, PyYAML, parc
R: data.table, dplyr, Matrix, future, future.apply, uwot, lisi, kBET, RANN, dpGMM, ggplot2
```

## External resources

The workflow requires several external resources, including the GRCh38 reference genome, RefSeq-based BED regions, Beagle genetic maps, a 1000 Genomes Project reference panel, ANNOVAR databases, and the CADD prescored database.

A complete list of required external resources, expected directory structure, download links, and preparation notes is provided in [**external resources guide**](docs/external_resources.md).

Paths to these resources should be specified in `config/local_config.yaml`, created from `config/example_config.yaml`.

## Input data

The workflow can start from either:

1. paired-end FASTQ files, followed by QC, trimming, alignment, duplicate marking, variant calling, and joint genotyping; or
2. existing BAM files, followed by optional liftover, duplicate-aware processing, variant calling, and joint genotyping.

Input filename conventions are controlled in `config/local_config.yaml`, including:

```yaml
fastq_r1_suffix: "_1.fastq"
fastq_r2_suffix: "_2.fastq"
deepvariant_bam_suffix: "_marked.bam"
```

The `deepvariant_bam` directory should point to the BAM files used directly by DeepVariant. For newly aligned data, this can be the duplicate-marked BAM directory. For lifted data, this can be the liftover output directory. The `deepvariant_bam_suffix` value must match the actual BAM suffix to avoid accidentally using unsorted, unmarked, or intermediate files.


## Repository structure

```text
config/                    # YAML configuration files
Data/                      # example inputs and auxiliary resources
Data_pre_processing/       # QC, trimming, alignment, duplicate marking, liftover
Variant_calling/           # DeepVariant and GLnexus
Variant_post_processing/   # genotype imputation, ANNOVAR/CADD annotation
Variant_to_gene/           # gene-level CWAF feature aggregation
Gene-level imputation/     # clustering, GMM thresholding, MNAR-aware feature imputation
docs/                      # installation and setup documentation
scripts/                   # helper scripts, including configuration parsing
```

## Workflow overview

### 1. Quality control, trimming, alignment, and liftover

Main scripts:

```text
Data_pre_processing/QC/trimm.sh
Data_pre_processing/Alignment/alignment.sh
Data_pre_processing/Liftover/liftover_bams.sh
```

This stage performs FASTQ quality control, optional adapter/quality trimming, BWA alignment to GRCh38, coordinate sorting, duplicate marking with Picard, and optional liftover of GRCh37/hg19-aligned BAM files to GRCh38/hg38.

### 2. Variant calling and joint genotyping

Main scripts:

```text
Variant_calling/Calling/run_deepvariant.sh
Variant_calling/Joint_genotyping/run_GLnexus.sh
```

DeepVariant is run per sample and produces VCF/GVCF files. GLnexus then merges the GVCF files and performs cohort-level joint genotyping. Both steps are restricted to the RefSeq exon/splice BED regions.

### 3. Genotype imputation

Main script:

```text
Variant_post_processing/1_genotype_imputation.sh
```

This stage normalizes variants, splits multiallelic sites, restricts records to SNVs/indels passing the configured allele-depth rule, conforms genotypes to the reference panel, phases haplotypes, and imputes missing genotypes with Beagle.

### 4. Variant annotation

Main script:

```text
Variant_post_processing/2_annotation.sh
```

Observed variants are annotated with ANNOVAR and CADD. The workflow keeps coding and splice-related records for downstream gene-level feature construction.

### 5. Gene-level feature generation

Main scripts:

```text
Variant_to_gene/gene_aggregation.sh
Variant_to_gene/cal_features_multi.py
```

Variants are aggregated by sample and gene. The main feature is CWAF, computed from allele fraction values weighted by CADD Phred-like scores. The procedure uses both the directly observed and genotype-imputed VCFs where applicable.

### 6. Feature analysis and gene-level imputation

Main scripts:

```text
Gene-level imputation/1_features_loading.R
Gene-level imputation/2_clustering.R
Gene-level imputation/3_GMM.R
Gene-level imputation/4_feature_imputation.R
```

This stage loads gene-level features, clusters samples with similar CWAF profiles, estimates detection-rate thresholds with Gaussian mixture modeling, defines MNAR entries, and performs masked cosine-similarity kNN imputation across technically comparable samples. Unlike for the remaining steps, the **user-adjustable parameters are kept at the beginning of these scripts, and can be modified if needed**.

## Reference panel construction
A detailed [**guide for the construction of reference panel**](docs/reference_panel.md)  is included in the repository. 

The benchmark EUR reference panel was constructed from 1000 Genomes Project high-coverage GRCh38 VCF files. The general procedure was:

1. select samples belonging to the EUR superpopulation from the 1kGP metadata;
2. process each chromosome separately;
3. keep biallelic SNPs and indels;
4. filter variants by minor allele count using `MAC>=5`;
5. remove `chr` prefixes from chromosome names to match the workflow reference convention;
6. index the resulting VCF files with tabix;
7. convert chromosome-specific VCF files to Beagle `bref3` format.

For chromosome X, only female EUR samples were used in the provided EUR panel preparation to avoid sex-ploidy inconsistencies during reference construction.

A helper script documenting this procedure is provided with the repository. The same logic can be adapted to construct a custom ancestry-matched panel by changing the sample-selection step in the metadata filter.

## Main outputs

The workflow creates output subdirectories under the configured `results_dir`, including intermediate genotype-imputation files, annotation outputs, feature matrices, and gene-level imputation results.

`results_dir` structure:

```text
Genotyping/        # GLnexus cohort-level outputs
Imputation/        # genotype-imputed VCFs and intermediate files
Annotation/        # ANNOVAR/CADD annotation outputs
Features/          # gene-level feature matrices
Results/           # downstream feature analysis and imputation results
Intermediate/      # temporary or step-specific intermediate files
```

## Notes and limitations

The workflow was designed for GRCh38/hg38 resources using chromosome names without the `chr` prefix. Mixing references, BED files, CADD files, and reference panels with inconsistent chromosome naming will cause downstream errors.

Large datasets such as FASTQ, BAM, VCF, CADD, and full reference-panel files are not distributed directly through GitHub. Users should place these files locally and provide paths in `config/local_config.yaml`.

The gene-level imputation strategy assumes that strong detection-rate differences across technical groups are unlikely to reflect true biology. For multi-ancestry cohorts, ancestry should be handled carefully because allele fraction-derived features may also reflect population structure.

