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

After installation and configuration are complete, users can run the lightweight example workflow described in the [example run guide](docs/example_run.md). The example uses four small FASTQ files reconstructed from selected regions of public SEQC2 WES BAM files. It is intended as a technical execution test, not as a biological benchmark. We highly reccomend the users to get familiar with example run, even if they do not plan to execute it, as it provides detailed explanation of workflow structure.

## Configuration

The workflow uses a YAML configuration file:

```text
config/example_config.yaml
```

Copy this file to a local, user-specific configuration file:

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

In the configuration file, `jobs` denotes externally parallelized processes, while `threads` denotes CPU threads used internally by a given tool. These values should be adapted to the available CPU cores, memory, storage throughput, and local scheduler limitations.

Some downstream analytical parameters are intentionally not stored in the global YAML file. Parameters controlling UMAP, clustering, detection-rate thresholding, and gene-level imputation are kept near the beginning of the corresponding R scripts. In particular, users should review and adjust the low- and high-detection-rate thresholds used for MNAR flagging in:

```text
Gene-level imputation/4_feature_imputation.R
```

The most important imputation parameters are:

```r
threshold_low_value <- 0.44
threshold_high_value <- 0.85
```

and were adjusted to match the [example run](docs/example_run.md). The low threshold defines when a gene is considered poorly detected in a given sample cluster, while the high threshold defines when the same gene is considered well detected in another sufficiently large cluster. These values control which missing gene-level CWAF values are treated as MNAR and therefore selected for imputation.


## Lightweight example run

We highly encourage the users to get familiar with our complete command-by-command  [**example run guide**](docs/example_run.md). 

The example workflow demonstrates the main processing steps from FASTQ files to gene-level feature generation. It is based on four small FASTQ files reconstructed from selected regions of BAM files and is intended only to verify installation, configuration, external resources, and basic workflow execution.

The public SEQC2 WES BAM files used to reconstruct the lightweight example FASTQ files originate from the benchmark dataset published by Zhao et al. (2021), which should be cited when using the example data.

Because the complete WES processing workflow is computationally intensive and produces large intermediate files, the example WES dataset is intentionally small. It is sufficient to validate the technical execution of the pipeline, but it is too small to meaningfully demonstrate cohort-level downstream analyses such as sample clustering, detection-rate modeling, MNAR masking, and gene-level imputation.

For this reason, the downstream gene-level imputation module is demonstrated in the example guide using a separate artificial feature-level dataset. This artificial dataset preserves the expected input structure and directory layout while providing enough samples to demonstrate clustering, detection-rate modeling, and imputation behavior.

## Workflow overview

### 1. Quality control, trimming, alignment, and liftover

Main scripts:

```text
Data_pre_processing/QC/trimm.sh
Data_pre_processing/Alignment/alignment.sh
Data_pre_processing/Liftover/liftover_bams.sh
```

This stage performs FASTQ quality control, optional adapter and quality trimming, BWA alignment to GRCh38, coordinate sorting, duplicate marking with Picard, and optional liftover of GRCh37/hg19-aligned BAM files to GRCh38/hg38.

User-adjustable parameters include FASTQ filename suffixes, Trimmomatic adapter and trimming settings, BWA alignment settings, the number of alignment jobs, thread counts, and the input/output directories.

### 2. Variant calling and joint genotyping

Main scripts:

```text
Variant_calling/Calling/run_deepvariant.sh
Variant_calling/Joint_genotyping/run_GLnexus.sh
```

DeepVariant is run per sample and produces VCF/GVCF files. GLnexus then merges the GVCF files and performs cohort-level joint genotyping. Both steps are restricted to the configured exon/splice BED regions.

User-adjustable parameters include the DeepVariant Docker image, BAM input directory, BAM suffix used for variant calling, number of shards, output directories, GLnexus input search root, and runtime resource settings.

### 3. Genotype imputation

Main script:

```text
Variant_post_processing/1_genotype_imputation.sh
```

This stage normalizes variants, splits multiallelic sites, restricts records to retained SNVs/indels, conforms genotypes to the reference panel, phases haplotypes, and imputes missing genotypes with Beagle.

Users should configure the Beagle and conform-gt paths, genetic maps, reference panel directory, chromosome naming convention, temporary directories, and parallelization settings. The workflow expects resources to use chromosome names compatible with the configured reference genome.

### 4. Variant annotation

Main script:

```text
Variant_post_processing/2_annotation.sh
```

Observed variants are annotated with ANNOVAR and CADD. The workflow keeps coding and splice-related records for downstream gene-level feature construction.

Users should configure ANNOVAR paths and databases, CADD script and prescored database paths, annotation output directories, and the number of parallel annotation jobs.

### 5. Gene-level feature generation

Main scripts:

```text
Variant_to_gene/gene_aggregation.sh
Variant_to_gene/cal_features_multi.py
```

Variants are aggregated by sample and gene. The main feature is CWAF, computed from allele fraction values weighted by CADD Phred-like scores. The procedure uses the directly observed CADD-scored VCF and, where available, the filtered genotype-imputed VCF.

If a filtered imputed VCF is not available for a chromosome, the aggregation step can calculate features from the observed CADD-scored VCF only. This behavior is useful for reduced example data and for cases where no variants are retained after genotype conformation or filtering.

### 6. Feature loading, clustering, detection-rate modeling, and gene-level imputation

Main scripts:

```text
Gene-level imputation/1_features_loading.R
Gene-level imputation/2_clustering.R
Gene-level imputation/3_GMM.R
Gene-level imputation/4_feature_imputation.R
```

This stage loads gene-level CWAF features, prepares sample-by-gene matrices, visualizes cohort structure with UMAP, clusters samples with similar CWAF profiles, models gene detection rates across clusters, defines MNAR candidate entries, and performs masked cosine-similarity kNN imputation.

The detection-rate modeling step is used to inspect the distribution of gene detection rates across sample clusters and to support the choice of low- and high-detection-rate thresholds. The thresholds used for MNAR flagging are defined manually at the beginning of `4_feature_imputation.R`, with default values as:
```r
threshold_low_value <- 0.44
threshold_high_value <- 0.85
```
selected to match the [**example run workflow**](docs/example_run.md).

A gene is considered potentially under-detected in a cluster when its detection rate in that cluster is below the low threshold and its detection rate in at least one other sufficiently large cluster is above the high threshold. Only missing CWAF values matching this MNAR pattern are selected for imputation.

The user should review these thresholds before applying the imputation step to a new cohort, because defaulte values were adjusted to artificially generated dataset. In reality, they may depend on cohort size, capture-kit composition, feature sparsity, and the observed detection-rate distribution.


## Main outputs

The workflow creates output subdirectories under the configured `results_dir`.

Typical output structure:

```text
Genotyping/              # GLnexus cohort-level outputs
Imputation/              # genotype-imputed VCFs and intermediate files
Annotation/              # ANNOVAR/CADD annotation outputs
Features/                # per-chromosome gene-level feature files
Gene_level_imputation/   # loaded feature matrices, UMAPs, clustering, GMM diagnostics, imputed matrices
Intermediate/            # temporary or step-specific intermediate files
```

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