# ![UCL-BLIC/rnaseq_variant_calling](docs/images/nfcore-rnaseq_logo.png)

[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)


### Introduction

**UCL-BLIC/rnaseq_variant_calling** is a bioinformatics analysis pipeline used for calling variants on RNA sequencing data

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)), aligns the reads 
 with [STAR](https://github.com/alexdobin/STAR) and does variant calling with [GATK](https://software.broadinstitute.org/gatk/) following the best practices [for RNAseq 
data](https://software.broadinstitute.org/gatk/documentation/article.php?id=3891)

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

### Documentation
The UCL-BLIC/rnaseq_variant_calling pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)
4. [Troubleshooting](docs/troubleshooting.md)

