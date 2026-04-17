# WESworkflow
A Practical Workflow for Correcting Kit-specific Effects in Whole-Exome Sequencing Data

---

## Requirements

* Bash, R (≥4.0), Python (for PARC)
* Tools: BWA, CrossMap (optional), samtools, Picard, DeepVariant, GLnexus, bcftools, Beagle, ANNOVAR
* Additional R packages listed in the R scripts.
* External resources:

  * Reference genome (GRCh38)
  * Reference panel (provided via Zenodo)
  * Genetic maps (for Beagle)
  * ANNOVAR databases
  * CADD prescored database

---

## Repository structure

```
Data/                      # example inputs and auxiliary resources
Data_pre_processing/       # QC + trimming + alignment / liftover
Variant_calling/           # DeepVariant + GLnexus
Variant_post_processing/   # genotype imputation with BEAGLE + ANNOVAR/CADD annotation
Variant_to_gene/           # gene-level variant aggregation
Gene-level imputation/     # clustering, GMM for threshold determination, gene-level imputation imputation
```

---

## Workflow

### 1. Data preprocessing

Scripts:

```
Data_pre_processing/QC/trimm.sh
Data_pre_processing/QC/fastqc/
Data_pre_processing/Alignment/alignment.sh
Data_pre_processing/Liftover/liftover_bams.sh
```

Steps:

* QC and trimming
* alignment to GRCh38
* duplicate marking
* optional liftover

---

### 2. Variant calling

Scripts:

```
Variant_calling/Calling/run_deepvariant.sh
Variant_calling/Joint_genotyping/run_GLnexus.sh
```

Steps:

* per-sample variant calling (DeepVariant)
* joint genotyping (GLnexus)

---

### 3. Genotype imputation

Script:

```
Variant_post_processing/1_genotype_imputation.sh
```

Steps:

* normalization and filtering
* genotype conformation
* imputation (Beagle)

---

### 4. Variant annotation

Script:

```
Variant_post_processing/2_annotation.sh
```

Steps:

* functional annotation (ANNOVAR)
* filtering to coding/splicing variants
* adding CADD scores

---

### 5. Gene-level feature generation

Scripts:

```
Variant_to_gene/gene_aggregation.sh
Variant_to_gene/cal_features_multi.py
```

Steps:

* aggregation of variants per gene
* computation of gene-level features

---

### 6. Feature analysis and imputation

Scripts:

```
Gene-level imputation/
├── 1_features_loading.R
├── 2_clustering.R
├── 3_GMM.R
├── 4_feature_imputation.R
```

Steps:

* feature loading
* clustering (PARC)
* detection rate modeling (GMM) for feature imputation thresholds
* MNAR-aware gene-level KNN imputation

---

## Data

* Example metadata:

```
Data/example/sample_path_map_example.tsv
```

* BED regions:

```
Data/bed/refGene_exons_splice5.nochr.bed
```

* LoF annotations:

```
Data/gnomad_lofs/
```

* Reference panel (EUR subset):

  * provided separately via Zenodo (see DOI)

---

## Notes

* All scripts use `/path/to/...` placeholders — paths must be adapted
* Large datasets (VCF, BAM, reference panels) are not included
* Workflow is modular — individual steps can be run independently
