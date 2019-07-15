# ![nfcore/rnaseq](docs/images/nfcore-rnaseq_logo.png)

[![Build Status](https://travis-ci.org/nf-core/rnaseq.svg?branch=master)](https://travis-ci.org/nf-core/rnaseq)
[![Nextflow](https://img.shields.io/badge/nextflow-%E2%89%A50.32.0-brightgreen.svg)](https://www.nextflow.io/)
[![DOI](https://zenodo.org/badge/127293091.svg)](https://zenodo.org/badge/latestdoi/127293091)
[![Gitter](https://img.shields.io/badge/gitter-%20join%20chat%20%E2%86%92-4fb99a.svg)](https://gitter.im/nf-core/Lobby)


### Introduction

**UCL-BLIC/rnaseq** is a bioinformatics analysis pipeline used for RNA sequencing data, modified to add kallisto.

The workflow processes raw data from FastQ inputs ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [Trim Galore!](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/)), aligns the reads 
([STAR](https://github.com/alexdobin/STAR) or [HiSAT2](https://ccb.jhu.edu/software/hisat2/index.shtml)), generates gene counts ([featureCounts](http://bioinf.wehi.edu.au/featureCounts/), 
[StringTie](https://ccb.jhu.edu/software/stringtie/)) as well as [kallisto](https://pachterlab.github.io/kallisto/) abundance files, and performs extensive quality-control on the results ([RSeQC](http://rseqc.sourceforge.net/), 
[dupRadar](https://bioconductor.org/packages/release/bioc/html/dupRadar.html), 
[Preseq](http://smithlabresearch.org/software/preseq/), [edgeR](https://bioconductor.org/packages/release/bioc/html/edgeR.html), [MultiQC](http://multiqc.info/)). See the [output documentation](docs/output.md) for more details of the results.

The pipeline is built using [Nextflow](https://www.nextflow.io), a bioinformatics workflow tool to run tasks across multiple compute infrastructures in a very portable manner.

### Documentation
The UCL-BLIC/rnaseq pipeline comes with documentation about the pipeline, found in the `docs/` directory:

1. [Installation](docs/installation.md)
2. [Running the pipeline](docs/usage.md)
3. [Output and how to interpret the results](docs/output.md)
4. [Troubleshooting](docs/troubleshooting.md)

### Credits
These scripts were originally written for use at the [National Genomics Infrastructure](https://portal.scilifelab.se/genomics/), part of [SciLifeLab](http://www.scilifelab.se/) in Stockholm, Sweden, by Phil Ewels ([@ewels](https://github.com/ewels)) and Rickard Hammarén ([@Hammarn](https://github.com/Hammarn)).

Many thanks to other who have helped out along the way too, including (but not limited to):
[@Galithil](https://github.com/Galithil),
[@pditommaso](https://github.com/pditommaso),
[@orzechoj](https://github.com/orzechoj),
[@apeltzer](https://github.com/apeltzer),
[@colindaven](https://github.com/colindaven).
