# UCL-BLIC/rnaseq_variant_calling Usage

## Running the pipeline
The typical command for running the pipeline is as follows:
```bash
module load blic-modules
module load nextflow_rnaseq_variant_calling

nextflow_rnaseq_variant_calling --reads '*_R{1,2}.fastq.gz' --genome hg38
```

This will launch the pipeline with the `legion` or `myriad` configuration profile, depending on where you submit the job from.

Note that the pipeline will create the following files in your working directory:

```bash
work            # Directory containing the nextflow working files
results         # Finished results (configurable, see below)
.nextflow_log   # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

## Main Arguments

### `-profile`
This parameter is NOT necessary as the shortcut `nextflow_rnaseq_variant_calling` takes care of selecting the appropiate configuration profile. But just for your information, profiles are used to give 
configuration presets for different compute environments.

* `legion`
    * A generic configuration profile to be used with the UCL cluster legion
* `myriad`
    * A generic configuration profile to be used with the UCL cluster myriad

### `--reads`
Use this to specify the location of your input FastQ files. For example:

```bash
--reads 'path/to/data/sample_*_{1,2}.fastq'
```

Please note the following requirements:

1. The path must be enclosed in quotes
2. The path must have at least one `*` wildcard character
3. When using the pipeline with paired end data, the path must use `{1,2}` notation to specify read pairs.

If left unspecified, a default pattern is used: `data/*{1,2}.fastq.gz`

### `--singleEnd`
By default, the pipeline expects paired-end data. If you have single-end data, you need to specify `--singleEnd` on the command line when you launch the pipeline. A normal glob pattern, enclosed in quotation marks, can then be used for `--reads`. For example:

```bash
--singleEnd --reads '*.fastq'
```

## Alignment tool
The pipeline uses [STAR](https://github.com/alexdobin/STAR) to align the raw FastQ reads to the reference genome. STAR is fast and common, but requires a lot of memory to run, typically around 38GB for the Human hg19 reference genome.

#### Minimum alignment lengths
When using STAR as the aligner then the pipeline sets the minimum alignment length to 15 base pairs. This filters out very short read alignments that can result from soft-clipping.

To customise this threshold, use the pipeline command line parameter `--min_aln_length [num]`, where `[num]` is the minimum required alignment length in base pairs. Mapping sensitivity can be improved by increasing this value, or made more flexible by decreasing it. The default value in STAR is `0`.

## Reference Genomes

The pipeline config files come bundled with paths to the following reference genome assemblies: hg19 and hg38. The pipeline has this aprameter set up as `false` by default, you need to specify it using the `--genome` 
flag.

* Human
  * `--genome hg19`
  * `--genome hg38`

### `--saveTrimmed`
By default, trimmed FastQ files will not be saved to the results directory. Specify this
flag (or set to true in your config file) to copy these files when complete.

### `--saveStarBAM`
As above, by default intermediate BAM files from the alignment will not be saved. The final BAM files created after the Recalibration step are always saved. Set to true to also copy out BAM files from STAR.

### `--saveDedupBAM`
As above, by default intermediate BAM files from the alignment will not be saved. The final BAM files created after the Recalibration step are always saved. Set to true to also copy out BAM files from MarkDuplicates (i.e., STAR + MarkDup).

## Adapter Trimming
If specific additional trimming is required (for example, from additional tags),
you can use any of the following command line parameters. These affect the command
used to launch TrimGalore!

### `--clip_r1 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads).

### `--clip_r2 [int]`
Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only).

### `--three_prime_clip_r1 [int]`
Instructs Trim Galore to remove bp from the 3' end of read 1 _AFTER_ adapter/quality trimming has been performed.

### `--three_prime_clip_r2 [int]`
Instructs Trim Galore to re move bp from the 3' end of read 2 _AFTER_ adapter/quality trimming has been performed.


## Library Prep Presets
Some command line options are available to automatically set parameters for common RNA-seq library preparation kits.

> Note that these presets override other command line arguments. So if you specify `--pico --clip_r1 0`, the `--clip_r1` bit will be ignored.

If you have a kit that you'd like a preset added for, please let us know!

### `--pico`
Sets trimming and standedness settings for the _SMARTer Stranded Total RNA-Seq Kit - Pico Input_ kit.

Equivalent to: `--forward_stranded` `--clip_r1 3` `--three_prime_clip_r2 3`

## Job Resources
### Automatic resubmission
Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the steps in the pipeline, if the job exits with an error code of `143` (exceeded requested resources) it will automatically resubmit with higher requests (2 x original, then 3 x original). If it still fails after three times then the pipeline is stopped.

###
## Other command line parameters
### `--outdir`
The output directory where the results will be saved.

### `-name`
Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.

This is used in the MultiQC report (if not default) and in the summary HTML / e-mail (always).

**NB:** Single hyphen (core Nextflow option)

### `-resume`
Specify this when restarting a pipeline. Nextflow will used cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously.

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

**NB:** Single hyphen (core Nextflow option)

