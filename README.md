# A practical workflow for correcting kit-specific effects in whole-exome sequencing data

**WESworkflow** is a modular pipeline for processing whole-exome sequencing (WES) data from raw reads or pre-aligned BAM files to gene-level variant feature matrices. The workflow performs quality control, trimming, alignment or liftover, duplicate marking, variant calling, joint genotyping, genotype imputation, functional annotation, CADD-weighted allele fraction (CWAF) aggregation, batch-effect assessment, and MNAR-aware gene-level imputation.

The workflow was developed for GRCh38/hg38 WES data and is intended for multi-source cohorts where exome capture kit differences may introduce systematic gene-level detection biases.

## Quick start

A detailed setup guide is provided in the [installation guide](docs/installation.md). A step-by-step lightweight execution example is provided in the [example run guide](docs/example_run.md).

In a typical setup, the user should:

```bash
git clone https://github.com/ZAEDPolSl/WESworkflow.git
cd WESworkflow
cp config/example_config.yaml config/local_config.yaml
```

Then edit `config/local_config.yaml` so that all tool paths, resource paths, input/output directories, and adjustable runtime parameters match the local environment.

Before running the full workflow on a cohort, run the lightweight installation check described in the [installation guide](docs/installation.md). This step verifies that dependencies, external binaries, reference files, and configuration paths are visible before launching computationally expensive jobs.



## Configuration

The workflow uses a YAML configuration file `config/example_config.yaml`. Copy this file to a local, user-specific configuration file:

```bash
cp config/example_config.yaml config/local_config.yaml
```


Then edit `config/local_config.yaml` so that all paths and runtime parameters match the local environment.

The most important settings that should be verified before running the workflow include:

```yaml
tools:
  deepvariant_docker_image: "google/deepvariant:1.8.0-gpu"
  conform_gt_jar: "/path/to/conform-gt.jar"
  beagle_jar: "/path/to/beagle.jar"
  annovar_dir: "/path/to/annovar"

resources:
  reference_fasta: "/path/to/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
  trimmomatic_adapters: "/path/to/TruSeq3-PE.fa"
  liftover_chain: "/path/to/GRCh37_to_GRCh38.chain.gz"
  refseq_bed: "Data/bed/refGene_exons_splice5.nochr.bed"
  genetic_maps_dir: "/path/to/beagle/genetic_maps"
  reference_panel_dir: "Data/reference_panel/EUR_nochr"
  cadd_prescored: "/path/to/CADD/whole_genome_SNVs.tsv.gz"

directories:
  raw_fastq: "/path/to/fastq"
  trimmed_fastq: "/path/to/trimmed_fastq_files"
  alignment_fastq: "/path/to/fastq_files"
  bam: "/path/to/aligned_bam_files"
  deepvariant_bam: "/path/to/deepvariant_input_bam_files"
  deepvariant_bam_suffix: "_marked.bam"
  deepvariant_output: "/path/to/vcf_files/snv"
  gvcf_search_root: "/path/to/general/datasets/parent/directory"
  results_dir: "/path/to/results"
  temporary_dir: "/path/to/tmp"
```

The `parameters` section contains user-adjustable runtime settings, including Trimmomatic parameters, FASTQ filename suffixes, BWA alignment settings, DeepVariant threads/GPU usage, genotype imputation settings, annotation parallelization, and gene-level feature aggregation parallelization.

The primary gene-level feature used in the workflow is the CADD-weighted average allele fraction (`CADD_weighted_avg_AF`, CWAF). An alternative cumulative metric (`CADD_weighted_cumulated_af`) is also available. The feature used for downstream analysis can be selected in the YAML configuration file:

```yaml
parameters:
  feature_column: "CADD_weighted_avg_AF"
```


In the configuration file, `jobs` denotes externally parallelized processes, while `threads` denotes CPU threads used internally by a given tool. These values should be adapted to the available CPU cores, memory, storage throughput, and local scheduler limitations.

Some downstream analytical parameters are intentionally **not stored** in the global YAML file but are defined at the beginning of the corresponding R scripts. While most parameters can be used with their default settings, particular attention should be paid to **PARC clustering and gene-level imputation parameters**, which should be adjusted to the characteristics of the analyzed cohort.
## Lightweight example run

We highly recommend reviewing the step-by-step [**example run guide**](docs/example_run.md) to become familiar with the workflow structure and execution.

The lightweight example uses four small FASTQ files reconstructed from public SEQC2 WES BAM files (Zhao et al., 2021) to verify installation, configuration, and basic pipeline execution. Due to its limited size, this dataset is not suitable for demonstrating cohort-level analyses.

Instead, downstream steps, including clustering, detection-rate modeling, and MNAR-aware gene-level imputation, are demonstrated using a separate artificial feature-level dataset with the expected input structure.

## Workflow overview

### 1. Quality control, trimming, alignment, and liftover

Main scripts:

```text
Data_pre_processing/QC/trimm.sh
Data_pre_processing/Alignment/alignment.sh
Data_pre_processing/Liftover/liftover_bams.sh
```

This stage performs FASTQ quality control, optional adapter and quality trimming, BWA alignment to GRCh38, coordinate sorting, duplicate marking with Picard, and optional liftover of GRCh37/hg19-aligned BAM files to GRCh38/hg38.


### 2. Variant calling and joint genotyping

Main scripts:

```text
Variant_calling/Calling/run_deepvariant.sh
Variant_calling/Joint_genotyping/run_GLnexus.sh
```

DeepVariant is run per sample and produces VCF/GVCF files. GLnexus then merges the GVCF files and performs cohort-level joint genotyping. Both steps are restricted to the configured exon/splice BED regions.



From this stage, all the downstream outputs will be generated in the directory specified under the `results_dir` parameter in the YAML file.

### 3. Genotype imputation

Main script:

```text
Variant_post_processing/1_genotype_imputation.sh
```

This stage normalizes variants, splits multiallelic sites, restricts records to retained SNVs/indels, conforms genotypes to the reference panel, phases haplotypes, and imputes missing genotypes with Beagle.



### 4. Variant annotation

Main script:

```text
Variant_post_processing/2_annotation.sh
```

Observed variants are annotated with ANNOVAR and CADD. The workflow keeps coding and splice-related SNV records for downstream gene-level feature construction.



### 5. Gene-level feature generation

Main scripts:

```text
Variant_to_gene/gene_aggregation.sh
Variant_to_gene/cal_features_multi.py
```

Variants are aggregated by sample and gene to generate CADD-weighted allele fraction (CWAF) features. The primary metric, `CADD_weighted_avg_AF`, represents the average allele fraction weighted by CADD Phred-like scores. An alternative metric, `CADD_weighted_cumulated_af`, additionally accounts for the cumulative contribution of variants within a gene.

Allele fractions are preferentially calculated from observed allelic depths (AD). When unavailable, confidently called homozygous-reference genotypes (GT = 0/0, DP ≥ 5) are assigned AF = 0; otherwise, imputed dosage (DS/2) is used when available.

The procedure combines information from the original CADD-scored VCF and, where available, the genotype-imputed VCF. If the filtered genotype-imputed VCF is unavailable for a chromosome, gene-level features are calculated using the original CADD-scored VCF only.
### 6. Feature loading, clustering, detection-rate modeling, and gene-level imputation

Main scripts:

```text
Gene-level imputation/1_features_loading.R
Gene-level imputation/2_clustering.R
Gene-level imputation/3_GMM.R
Gene-level imputation/4_feature_imputation.R
```

This stage loads the selected gene-level feature, constructs sample-by-gene matrices, and visualizes cohort structure using UMAP. PARC clustering identifies groups of samples with similar feature profiles reflecting the technical structure, while Gaussian mixture modeling (GMM) of gene detection rates provides candidate thresholds for identifying missing-not-at-random (MNAR) values. These values are subsequently imputed using masked cosine-similarity kNN.

#### 6.1 Parameter selection

Particular attention should be paid to the following parameters:

a) **PARC clustering (`2_clustering.R`):** `knn` and `resolution` should be adjusted based on the observed sample cluster structure.
  - *Users can run `2_clustering.R` multiple times with different `knn` and `resolution` settings, visually comparing the resulting PARC cluster assignments on UMAP. Each run saves a separate cluster mapping file in `Results/Clustering/`, with the corresponding parameters included in its filename. Once the optimal configuration has been selected, the corresponding file should replace `Results/sample_kit_cluster_map.tsv` before proceeding to GMM modeling and gene-level imputation.*




b) **Gene-level imputation (`4_feature_imputation.R`):** `threshold_low_value` and `threshold_high_value` should be selected based on the GMM results from `3_GMM.R`.


The current default settings were adjusted for the [example run](docs/example_run.md) and should be reviewed and adapted to the characteristics and results of each analyzed cohort.

## Main outputs

The workflow creates output subdirectories under the configured `results_dir`.

Typical output structure:

```text
<results_dir>/
    Genotyping/                 # GLnexus cohort-level outputs
    Imputation/                 # Genotype-imputed VCFs and intermediate files
    Annotation/                 # ANNOVAR/CADD annotation outputs
    Features/                   # Per-chromosome gene-level feature files
    Gene_level_imputation/      # Gene-level analysis for each selected feature
        <feature_column>/
            Results/            # Feature matrices, clustering, GMM, and gene-level imputation results
                Clustering/     # Cluster assignments for different PARC configurations
            Figures/            # UMAP visualizations and clustering plots
    Intermediate/               # Temporary or step-specific intermediate files
```

Gene-level analysis outputs are organized separately for each selected feature (`feature_column`). Some output filenames include the corresponding analysis parameters to distinguish different configurations.

For the lightweight example run, outputs are written under:

```text
Data/example/output/results/
```

A detailed list of expected files for each step is provided in the [example run guide](docs/example_run.md).

## Notes and limitations

The workflow was designed for GRCh38/hg38 resources using chromosome names without the `chr` prefix. Mixing references, BED files, CADD files, and reference panels with inconsistent chromosome naming will cause downstream errors.

Large datasets such as FASTQ, BAM, VCF, CADD, and full reference-panel files are not distributed directly through GitHub. Users should place these files locally and provide paths in `config/local_config.yaml`.

The lightweight example run is intended to verify workflow execution and configuration. It should not be interpreted as a biological WES analysis or as a benchmark of batch-effect correction.

The artificial downstream dataset used in the example run is provided only to demonstrate expected input structure, clustering, detection-rate modeling, MNAR masking, and gene-level imputation behavior. It is not intended to represent a biologically meaningful WES cohort.

The gene-level imputation strategy assumes that strong detection-rate differences across technical groups are unlikely to reflect true biology. For multi-ancestry cohorts, ancestry should be handled carefully because allele fraction-derived features may also reflect population structure.

## References

Zhao, Y., Fang, L.T., Shen, T.W., et al. Whole genome and exome sequencing reference datasets from a multi-center and cross-platform benchmark study. *Scientific Data* 8, 296 (2021). https://doi.org/10.1038/s41597-021-01077-5