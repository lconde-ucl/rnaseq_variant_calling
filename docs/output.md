# UCL-BLIC/rnaseq_variant_calling Output

UCL-BLIC/rnaseq_variant_calling is a pipeline to call variants on RNAseq data. This document describes the output produced by the pipeline.

Most of the plots are taken from the MultiQC report, which summarises results at the end of the pipeline.

## Pipeline overview
The pipeline is built using [Nextflow](https://www.nextflow.io/)
and processes data using the following steps:

* [FastQC](#fastqc) - read quality control
* [TrimGalore](#trimgalore) - adapter trimming
* [STAR](#star) - alignment
* [Picard](#picard) - mark PCR duplicates
* [GATK](#gatk) - Genome Analysis Toolkit
	- [Picard AddRGs](#addRGs) - add read groups
	- [Split & Trim](#split&trim) - split and trim
	- [BAM recalibration](#bam-recalibration) - recalibrate BAM file
	- [HaplotypeCaller](#haplotypecaller) - call variants in GVCF mode
	- [GenotypeGVCFs](#genotypegvcfs) - genotype generate GVCFs
	- [FilterVariants](#filtervariants) - hard filter variants
* [Annovar](#annovar) - annotate variants
* [MultiQC](#multiqc) - aggregate report, describing results of the whole pipeline

## FastQC
[FastQC](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/) gives general quality metrics about your reads. It provides information about the quality score distribution across your reads, the per base sequence content (%T/A/G/C). You get information about adapter contamination and other overrepresented sequences.

For further reading and documentation see the [FastQC help](http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/).

> **NB:** The FastQC plots displayed in the MultiQC report shows _untrimmed_ reads. They may contain adapter sequence and potentially regions with low quality. To see how your reads look after trimming, look at the FastQC reports in the `trim_galore` directory.

**Output directory: `results/fastqc`**

* `sample_fastqc.html`
  * FastQC report, containing quality metrics for your untrimmed raw fastq files
* `zips/sample_fastqc.zip`
  * zip file containing the FastQC report, tab-delimited data file and plot images

## TrimGalore
The pipeline uses [TrimGalore](http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) for removal of adapter contamination and trimming of low quality regions. TrimGalore uses [Cutadapt](https://github.com/marcelm/cutadapt) for adapter trimming and runs FastQC after it finishes.

MultiQC reports the percentage of bases removed by TrimGalore in the _General Statistics_ table, along with a line plot showing where reads were trimmed.

**Output directory: `results/trim_galore`**

Contains FastQ files with quality and adapter trimmed reads for each sample, along with a log file describing the trimming.

* `sample_val_1.fq.gz`, `sample_val_2.fq.gz`
  * Trimmed FastQ data, reads 1 and 2.
  * NB: Only saved if `--saveTrimmed` has been specified.
* `logs/sample_val_1.fq.gz_trimming_report.txt`
  * Trimming report (describes which parameters that were used)
* `FastQC/sample_val_1_fastqc.zip`
  * FastQC report for trimmed reads

Single-end data will have slightly different file names and only one FastQ file per sample.

## STAR
STAR is a read aligner designed for RNA sequencing.  STAR stands for Spliced Transcripts Alignment to a Reference, it produces results comparable to TopHat (the aligned previously used by NGI for RNA alignments) but is much faster.

The STAR section of the MultiQC report shows a bar plot with alignment rates: good samples should have most reads as _Uniquely mapped_ and few _Unmapped_ reads.

![STAR](images/star_alignment_plot.png)

**Output directory: `results/STAR`**

* `Sample_Aligned.sortedByCoord.out.bam`
  * The aligned BAM file
* `Sample_Log.final.out`
  * The STAR alignment report, contains mapping results summary
* `Sample_Log.out` and `Sample_Log.progress.out`
  * STAR log files, containing a lot of detailed information about the run. Typically only useful for debugging purposes.
* `Sample_SJ.out.tab`
  * Filtered splice junctions detected in the mapping

## Picard
[Picard MarkDuplicates](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates) locates and tags duplicate reads in a BAM or SAM file, where duplicate reads are defined as originating from a single fragment of DNA. 
Duplicates can arise during sample preparation e.g. library construction using PCR.

The MarkDuplicates tool works by comparing sequences in the 5 prime positions of both reads and read-pairs in a SAM/BAM file. After duplicate reads are collected, the tool differentiates the primary and duplicate reads using an algorithm that 
ranks reads by the sums of their base-quality scores (default method). For more information visit the [picard docs](http://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates).

**Output directory: `results/markDuplicates`**


## GATK

* [GATK](#gatk) - Genome Analysis Toolkit
	- [AddRGs](#addRGs) - add read groups
	- [Split & Trim](#split&trim) - split and trim
	- [BAM recalibration](#bam-recalibration) - recalibrate BAM file
	- [HaplotypeCaller](#haplotypecaller) - call variants in GVCF mode
	- [GenotypeGVCFs](#genotypegvcfs) - genotype generate GVCFs
	- [FilterVariants](#filtervariants) - hard filter variants

The pipeline uses the suggested [GATK](https://software.broadinstitute.org/gatk/) pipeline for [variant calling on RNAseq data] (https://software.broadinstitute.org/gatk/documentation/article.php?id=3891), and follows the following steps:

* picard AddOrReplaceReadGroups: adding read group information
* SplitNCigarReads: Next, we use a GATK tool called SplitNCigarReads developed specially for RNAseq, which splits reads into exon segments (getting rid of Ns but maintaining grouping information) and hard-clip any sequences overhanging into the 
intronic regions.
* BaseRecalibrator: The GATK [BaseRecalibrator](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_bqsr_BaseRecalibrator.php) detects systematic errors in 
base quality scores.
* HaplotypeCaller: The GATK [HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php) calls germline SNPs, insertions and deletions via 
local re-assembly of haplotypes.
* GenotypeGVCFs: The GATK [GenotypeGVCFs]() performs joint genotyping on gVCF files produced by the GATK 
[HaplotypeCaller](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php).
* VariantFiltration: Filter the resulting callset using hard filters

**Output directory: `results/variants`**

## Annovar

[annovar](http://annovar.openbioinformatics.org/en/latest) is a variant annotation and effect prediction tool. It annotates and predicts the effects of genetic variants (such as amino acid changes). 

**Output directory: `results/variants`**

## MultiQC
[MultiQC](http://multiqc.info) is a visualisation tool that generates a single HTML report summarising all samples in your project. Most of the pipeline QC results are visualised in the report and further statistics are available in within the 
report data directory.

**Output directory: `results/MultiQC`**

* `Project_multiqc_report.html`
  * MultiQC report - a standalone HTML file that can be viewed in your web browser
* `Project_multiqc_data/`
  * Directory containing parsed statistics from the different tools used in the pipeline

For more information about how to use MultiQC reports, see http://multiqc.info
