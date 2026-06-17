# External resources

This workflow requires several external resources that are not fully stored in the repository because of file size, licensing restrictions, or because they should be downloaded directly from the original provider. This document lists all resources required to run the workflow and describes the expected directory structure.

All paths to external resources should be configured in `config/example_config.yaml`.

## Genome build

The current workflow is configured for GRCh38 / hg38 resources. All reference files, annotation databases, BED files, genetic maps, CADD files, and reference panels should correspond to GRCh38 / hg38 unless explicitly stated otherwise.

The workflow assumes chromosome names without the `chr` prefix. If downloaded resources use chromosome names with the `chr` prefix, they should be converted before use.

## Expected directory structure

The following directory structure is recommended:

```
WESworkflow/
├── config/
│   └── example_config.yaml
├── Data/
│   ├── bed/
│   │   └── refGene_exons_splice5.nochr.bed
│   ├── reference_panel/
│   │   └── EUR_nochr/
│   ├── gnomad_lofs/
│   └── example/
└── resources/
    ├── GRCh38/
    ├── CADD/
    ├── beagle/
    │   └── genetic_maps/
    └── annovar/
        └── humandb/
```

The `resources/` directory is only a suggested location. Users may store large external resources elsewhere, but the corresponding paths must be updated in `config/example_config.yaml`.

## Resources included in the repository

The following lightweight resources are provided directly in the repository:

| Resource             | Location                                   | Notes                                                                                                                   |
| -------------------- | ------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------- |
| RefSeq exon BED file | `Data/bed/refGene_exons_splice5.nochr.bed` | BED file with RefSeq exon coordinates and 5 bp splice-region flanks, prepared for GRCh38/hg38 without the `chr` prefix. |
| Example metadata     | `Data/example/sample_path_map_example.tsv` | Example input table showing the expected metadata/sample path format.                                                   |
| LoF annotation files | `Data/gnomad_lofs/`                        | Auxiliary LoF annotation files used during feature calculation, if applicable.                                          |

Large resources such as reference genomes, CADD files, ANNOVAR databases, Beagle maps, and genotype imputation reference panels are not stored directly in the repository.

## Required external resources

| Resource                    | Required version/build                         | Config key                       | Expected file or directory                                                  | Source                          |
| --------------------------- | ---------------------------------------------- | -------------------------------- | --------------------------------------------------------------------------- | ------------------------------- |
| Reference genome FASTA      | GRCh38 / hg38                                  | `resources.reference_fasta`      | FASTA file with `.fai` and `.dict` indexes                                  | Ensembl GRCh38 primary assembly |
| Liftover chain file         | GRCh37/hg19 to GRCh38/hg38                     | `resources.liftover_chain`       | `.chain.gz` file                                                            | Ensembl assembly mapping        |
| RefSeq BED file             | hg38, no `chr` prefix                          | `resources.refseq_bed`           | `refGene_exons_splice5.nochr.bed`                                           | Provided in repository          |
| Beagle genetic maps         | GRCh38-compatible                              | `resources.genetic_maps_dir`     | Directory with chromosome-specific genetic map files                        | Beagle website                  |
| Genotype imputation reference panel             | 1000 Genomes Project, GRCh38, no `chr` prefix  | `resources.reference_panel_dir`  | Directory with chromosome-specific VCF/BCF files and indexes                | Provided via Zenodo                          |
| ANNOVAR                     | Current version compatible with hg38 databases | `tools.annovar_dir`              | ANNOVAR installation directory containing `table_annovar.pl` and `humandb/` | ANNOVAR website                 |
| ANNOVAR databases           | hg38                                           | `tools.annovar_dir/humandb`      | Required annotation databases inside `humandb/`                             | ANNOVAR website                 |
| CADD prescored SNV database | CADD v1.7 or compatible hg38 release           | `resources.cadd_prescored`       | `.tsv.gz` file with corresponding `.tbi` index                              | CADD download page              |
| DeepVariant Docker image    | v1.8.0                                         | `tools.deepvariant_docker_image` | Docker image name                                                           | DeepVariant GitHub repository             |
| Beagle JAR                  | Beagle v5.5 or compatible                      | `tools.beagle_jar`               | `.jar` file                                                                 | Beagle website                  |
| conform-gt JAR              | Compatible with Beagle preprocessing           | `tools.conform_gt_jar`           | `.jar` file                                                                 | Beagle website                  |


## Reference genome
The workflow is configured for the GRCh38 / hg38 genome build. The recommended reference is the Ensembl GRCh38.p14 primary assembly FASTA, available from the [Ensembl archive](http://feb2023.archive.ensembl.org/index.html). Because the workflow is restricted to canonical chromosomal regions and does not rely on alternative contigs or patch-specific sequences, other GRCh38 primary assembly releases may also be compatible, provided that chromosome naming and coordinate conventions are consistent with the remaining resources.

The reference FASTA should be indexed before running the workflow. The following files are expected next to the FASTA file:

| File            | Purpose                                                           |
| --------------- | ----------------------------------------------------------------- |
| `.fai`          | FASTA index generated by Samtools.                                |
| `.dict`         | Sequence dictionary required by Picard and some downstream tools. |
| BWA index files | Required for read alignment with BWA.                             |

Example expected files:

```
Homo_sapiens.GRCh38.dna.primary_assembly.fa
Homo_sapiens.GRCh38.dna.primary_assembly.fa.fai
Homo_sapiens.GRCh38.dna.primary_assembly.dict
Homo_sapiens.GRCh38.dna.primary_assembly.fa.amb
Homo_sapiens.GRCh38.dna.primary_assembly.fa.ann
Homo_sapiens.GRCh38.dna.primary_assembly.fa.bwt
Homo_sapiens.GRCh38.dna.primary_assembly.fa.pac
Homo_sapiens.GRCh38.dna.primary_assembly.fa.sa
```

Example commands:

```
samtools faidx Homo_sapiens.GRCh38.dna.primary_assembly.fa
bwa index Homo_sapiens.GRCh38.dna.primary_assembly.fa
picard CreateSequenceDictionary R=Homo_sapiens.GRCh38.dna.primary_assembly.fa O=Homo_sapiens.GRCh38.dna.primary_assembly.dict
```

## Liftover chain file

If BAM files aligned to GRCh37/hg19 are used as input, the workflow can lift coordinates over to GRCh38/hg38 using CrossMap. This requires a `GRCh37_to_GRCh38.chain.gz` chain file, which can be obtained from [Ensembl FTP site](https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/).

The chain file path should be set in:

```
resources.liftover_chain
```

If all input data are already aligned to GRCh38/hg38, this resource is not required.

## RefSeq BED file

Variant calling and joint genotyping are restricted to RefSeq-derived exonic regions with 5 bp splice-region flanks. The processed BED file is included in the repository:

```
Data/bed/refGene_exons_splice5.nochr.bed
```

The file is prepared for GRCh38/hg38 and uses chromosome names without the `chr` prefix. Script used for generating the regions and download link for corresponding RefSeq database can be found at `Variant_calling/Regions/get_target_regions.sh`

The corresponding configuration key is:

```
resources.refseq_bed
```

## Beagle genetic maps

Genotype imputation with Beagle requires chromosome-specific genetic maps which can be obtained from [Beagle website](https://faculty.washington.edu/browning/beagle/beagle.html). The directory containing these files should be set in:

```
resources.genetic_maps_dir
```

The files should be compatible with the genome build used by the reference panel and query VCF files. 

## Reference panel

The genotype imputation reference panel used in the benchmark analysis was constructed from the [1000 Genomes Project high-coverage GRCh38](https://www.internationalgenome.org/data-portal/data-collection/1000genomes_30x) VCF files and restricted to individuals of European ancestry. However, users can adjust the workflow to any ancestry available within the reference cohort by following the `docs/reference_panel.md` tutorial.

Due to file size, the prepared European reference panel used in this workflow is provided separately through Zenodo repository under the accession number [19626919](https://doi.org/10.5281/zenodo.19626919).

The expected location is:

```
Data/reference_panel/EUR_nochr/
```

The corresponding configuration key is:

```
resources.reference_panel_dir
```

The reference panel should contain chromosome-specific files and their indexes. File names should be consistent with the imputation scripts.

Example expected structure:

```
Data/reference_panel/EUR_nochr/
├── 1kGP.1.EUR.nochr.vcf.vcf.gz
├── 1kGP.1.EUR.nochr.vcf.vcf.gz.tbi
├── 1kGP.2.EUR.nochr.vcf.vcf.gz
├── 1kGP.2.EUR.nochr.vcf.vcf.gz.tbi
├── ...
└── BREF
    ├── 1kGP.1.EUR.nochr.bref3
    ├── 1kGP.2.EUR.nochr.bref3
    ├── ...

```

If the scripts expect chromosome names without the `chr` prefix, the panel files and VCF records should be prepared accordingly.

## ANNOVAR annotation databases

Variant annotation is performed with ANNOVAR. The ANNOVAR installation directory should contain `table_annovar.pl` and the `humandb/` directory. All databases should be downloaded for the hg38 genome build. Both the tool and databases are available at [ANNOVAR download](https://annovar.openbioinformatics.org/en/latest/).

The corresponding configuration key is:

```
tools.annovar_dir
```

The following ANNOVAR databases were used in the benchmark analysis:

| Database      | Genome build | Version    |
| ------------- | ------------ | ---------- |
| RefGene       | hg38         | 2020-08-17 |
| ExAC          | hg38         | v0.3       |
| gnomAD exomes | hg38         | v4.1       |
| 1000 Genomes  | hg38         | 2015-08    |
| dbSNP         | hg38         | build 150  |
| ClinVar       | hg38         | 2024-12-1  |
| dbNSFP        | hg38         | v4.7a      |

Expected structure:

```
annovar/
├── table_annovar.pl
└── humandb/
    ├── hg38_refGene.txt
    ├── hg38_refGeneMrna.fa
    ├── hg38_exac03.txt
    ├── hg38_gnomad41_exome.txt
    ├── hg38_ALL.sites.2015_08.txt
    ├── hg38_avsnp150.txt
    ├── hg38_clinvar_20241215.txt
    └── hg38_dbnsfp47a.txt
```


## CADD prescored database

The workflow uses the CADD v1.7 prescored SNV database during variant annotation and gene-level feature calculation, which can be obtained directly for the [CADD downloads website](https://cadd.bihealth.org/download)

The corresponding configuration key is:

```
resources.cadd_prescored
```

Expected files:

```
whole_genome_SNVs.tsv.gz
whole_genome_SNVs.tsv.gz.tbi
```

The `.tbi` index must be present next to the compressed CADD file.

## Resource checks before running the workflow

Before running the workflow, users should verify that all paths in `config/example_config.yaml` point to existing files or directories.

At minimum, check the following:

| Config key                       | Expected type                             |
| -------------------------------- | ----------------------------------------- |
| `resources.reference_fasta`      | file                                      |
| `resources.liftover_chain`       | file, only if liftover is used            |
| `resources.refseq_bed`           | file                                      |
| `resources.genetic_maps_dir`     | directory                                 |
| `resources.reference_panel_dir`  | directory                                 |
| `resources.cadd_prescored`       | file                                      |
| `tools.annovar_dir`              | directory                                 |
| `tools.beagle_jar`               | file                                      |
| `tools.conform_gt_jar`           | file                                      |

## Notes on reproducibility

To reproduce the benchmark analysis, users should use the same genome build, RefSeq BED file, CADD release, ANNOVAR database versions, and reference panel as described in the manuscript.

Using different genome builds, chromosome naming conventions, CADD releases, or reference panels may change the set of retained variants and therefore the final gene-level feature matrix.
