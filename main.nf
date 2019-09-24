/*
===============================================================
 UCL-BLIC/rnaseq_variant_calling
===============================================================
 RNA-Seq Variant Calling Pipeline. Started June 2019.
 #### Homepage / Documentation
 https://github.com/UCL-BLIC/rnaseq_Variant_calling
---------------------------------------------------------------
*/

def helpMessage() {
    log.info """
    =======================================================
                                              ,--./,-.
              ___     __   __   __   ___     /,-._.--~\'
        |\\ | |__  __ /  ` /  \\ |__) |__         }  {
        | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                              `._,._,\'
     RNA-Seq Variant Calling 
    =======================================================
    Usage:
    The typical command for running the pipeline is as follows:
    nextflow_rnaseq_variant_calling --reads '*_R{1,2}.fastq.gz' --genome GRCh37
    Mandatory arguments:
      --reads                       Path to input data (must be surrounded with quotes)
      --genome                      Name of iGenomes reference
      -profile                      Configuration profile to use. myriad / legion
    Options:
      --singleEnd                   Specifies that the input is single end reads
    Intermediate files
      --saveTrimmed                 Save trimmed FastQ file intermediates
      --saveStarBAM                 Save the BAM files from the Alignment step (after STAR) - not done by default
      --saveDedupBAM                Save dedup BAM (after STAR and markDup) - not saved by default
    Trimming options
      --clip_r1 [int]               Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2 [int]               Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1 [int]   Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2 [int]   Instructs Trim Galore to re move bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
    Other options:
      --outdir                      The output directory where the results will be saved
      -w/--work-dir                 The temporary directory where intermediate data will be saved
      -name                         Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic.
    """.stripIndent()
}

/*
 * SET UP CONFIGURATION VARIABLES
 */


// Show help emssage
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

// Check if genome exists in the config file
if (params.genomes && params.genome && !params.genomes.containsKey(params.genome)) {
    exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
  }

// Configurable variables
params.name = false
params.genome = false
params.star_index = params.genome ? params.genomes[ params.genome ].star ?: false : false
params.fasta = params.genome ? params.genomes[ params.genome ].fasta ?: false : false
params.gtf = params.genome ? params.genomes[ params.genome ].gtf ?: false : false
params.dbsnp = params.genome ? params.genomes[ params.genome ].dbsnp ?: false : false
params.annovar_anno = params.genome ? params.genomes[ params.genome ].annovar_anno ?: false : false


multiqc_config = file("$baseDir/assets/multiqc_config.yaml")
output_docs = file("$baseDir/docs/output.md")
wherearemyfiles = file("$baseDir/assets/where_are_my_files.txt")

// Define regular variables so that they can be overwritten
clip_r1 = params.clip_r1
clip_r2 = params.clip_r2
three_prime_clip_r1 = params.three_prime_clip_r1
three_prime_clip_r2 = params.three_prime_clip_r2


// Has the run name been specified by the user?
//  this has the bonus effect of catching both -name and --name
custom_runName = params.name
if( !(workflow.runName ==~ /[a-z]+_[a-z]+/) ){
  custom_runName = workflow.runName
}


/*
 * Create a channel for input read files
 */
if(params.readPaths){
    if(params.singleEnd){
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_trimgalore; }
    } else {
        Channel
            .from(params.readPaths)
            .map { row -> [ row[0], [file(row[1][0]), file(row[1][1])]] }
            .ifEmpty { exit 1, "params.readPaths was empty - no input files supplied" }
            .into { raw_reads_fastqc; raw_reads_trimgalore; }
    }
} else {
    Channel
        .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
        .ifEmpty { exit 1, "Cannot find any reads matching: ${params.reads}\nNB: Path needs to be enclosed in quotes!\nNB: Path requires at least one * wildcard!\nIf this is single-end data, please specify --singleEnd on the command line." }
        .into { raw_reads_fastqc; raw_reads_trimgalore; }
}


log.info """=======================================================
                                          ,--./,-.
          ___     __   __   __   ___     /,-._.--~\'
    |\\ | |__  __ /  ` /  \\ |__) |__         }  {
    | \\| |       \\__, \\__/ |  \\ |___     \\`-._,-`-,
                                          `._,._,\'
    RNA-Seq Variant Calling
======================================================="""
def summary = [:]
summary['Run Name']     = custom_runName ?: workflow.runName
summary['Reads']        = params.reads
summary['Data Type']    = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Genome']       = params.genome
summary['Trim R1'] = clip_r1
summary['Trim R2'] = clip_r2
summary["Trim 3' R1"] = three_prime_clip_r1
summary["Trim 3' R2"] = three_prime_clip_r2
summary['Aligner'] = "STAR"
summary['STAR Index']   = params.star_index
summary['Fasta Ref']    = params.fasta
summary['GTF Annotation']  = params.gtf
summary['Save Trimmed']   = params.saveTrimmed ? 'Yes' : 'No'
summary['Save Aligned BAM'] = params.saveStarBAM ? 'Yes' : 'No'
summary['Save markDuplicated BAM'] = params.saveDedupBAM ? 'Yes' : 'No'
summary['Max Memory']     = params.max_memory
summary['Max CPUs']       = params.max_cpus
summary['Max Time']       = params.max_time
summary['Output dir']     = params.outdir
summary['Working dir']    = workflow.workDir
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current home']   = "$HOME"
summary['Current user']   = "$USER"
summary['Current path']   = "$PWD"
summary['Script dir']     = workflow.projectDir
summary['Config Profile'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(15)}: $v" }.join("\n")
log.info "========================================="




/*
 * STEP 1 - FastQC
 */
process fastqc {
    tag "$name"
    publishDir "${params.outdir}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(reads) from raw_reads_fastqc

    output:
    file "*_fastqc.{zip,html}" into fastqc_results

    script:
    """
    fastqc -q $reads
    """
}


/*
 * STEP 2 - Trim Galore!
 */
process trim_galore {
    tag "$name"
    publishDir "${params.outdir}/trim_galore", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
            else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
            else if (!params.saveTrimmed && filename == "where_are_my_files.txt") filename
            else if (params.saveTrimmed && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    set val(name), file(reads) from raw_reads_trimgalore
    file wherearemyfiles

    output:
    file "*fq.gz" into trimmed_reads, trimmed_reads_kallisto
    file "*trimming_report.txt" into trimgalore_results
    file "*_fastqc.{zip,html}" into trimgalore_fastqc_reports
    file "where_are_my_files.txt"


    script:
    c_r1 = clip_r1 > 0 ? "--clip_r1 ${clip_r1}" : ''
    c_r2 = clip_r2 > 0 ? "--clip_r2 ${clip_r2}" : ''
    tpc_r1 = three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${three_prime_clip_r1}" : ''
    tpc_r2 = three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${three_prime_clip_r2}" : ''
    if (params.singleEnd) {
        """
        trim_galore --fastqc --gzip $c_r1 $tpc_r1 $reads
        """
    } else {
        """
        trim_galore --paired --fastqc --gzip $c_r1 $c_r2 $tpc_r1 $tpc_r2 $reads
        """
    }
}


/*
 * STEP 3 - align with STAR
 */
// Function that checks the alignment rate of the STAR output
// and returns true if the alignment passed and otherwise false
skipped_poor_alignment = []
def check_log(logs) {
    def percent_aligned = 0;
    logs.eachLine { line ->
        if ((matcher = line =~ /Uniquely mapped reads %\s*\|\s*([\d\.]+)%/)) {
            percent_aligned = matcher[0][1]
        }
    }
    logname = logs.getBaseName() - 'Log.final'
    if(percent_aligned.toFloat() <= '5'.toFloat() ){
        log.info "#################### VERY POOR ALIGNMENT RATE! IGNORING FOR FURTHER DOWNSTREAM ANALYSIS! ($logname)    >> ${percent_aligned}% <<"
        skipped_poor_alignment << logname
        return false
    } else {
        log.info "          Passed alignment > star ($logname)   >> ${percent_aligned}% <<"
        return true
    }
}

process star {
    tag "$prefix"
    publishDir "${params.outdir}/STAR", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf(".bam") == -1) "logs/$filename"
            else if (!params.saveStarBAM && filename == "where_are_my_files.txt") filename
            else if (params.saveStarBAM && filename != "where_are_my_files.txt") filename
            else null
        }
    
    input:
    file reads from trimmed_reads
    file wherearemyfiles

    output:
    set file("*Log.final.out"), file ('*.bam') into star_aligned
    file "*.out" into alignment_logs
    file "*SJ.out.tab"
    file "where_are_my_files.txt"

    script:
    prefix = reads[0].toString() - ~/(_R1)?(_001)?(_trimmed)?(_val_1)?(\.fq)?(\.fastq)?(\.gz)?$/
    def avail_mem = task.memory ? "--limitBAMsortRAM ${task.memory.toBytes() - 100000000}" : ''
    """
    STAR --genomeDir ${params.star_index} \\
        --sjdbGTFfile ${params.gtf} \\
        --readFilesIn $reads  \\
        --runThreadN ${task.cpus} \\
        --twopassMode Basic \\
        --outWigType bedGraph \\
        --outSAMtype BAM SortedByCoordinate $avail_mem \\
        --readFilesCommand zcat \\
        --runDirPerm All_RWX \\
        --outFileNamePrefix $prefix \\
    """
}
// Filter removes all 'aligned' channels that fail the check
star_aligned
    .filter { logs, bams -> check_log(logs) }
    .flatMap {  logs, bams -> bams }
.into { bam_markduplicates }




/*
 * STEP 6 Mark duplicates
 */
process markDuplicates {

    tag "${bam.baseName - 'Aligned.sortedByCoord.out'}"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy',
        saveAs: {filename ->
            if (filename.indexOf("_metrics.txt") > 0) "metrics/$filename"
            else if (!params.saveDedupBAM && filename == "where_are_my_files.txt") filename
            else if (params.saveDedupBAM && filename != "where_are_my_files.txt") filename
            else null
        }

    input:
    file bam from bam_markduplicates
    file wherearemyfiles

    output:
    set val("${sample}"), file("${sample}.markDups.bam"), file("${sample}.markDups.bai") into bam_md
    file "${sample}.markDups_metrics.txt" into picard_resuls
    file "where_are_my_files.txt"

    script:
    sample = bam.baseName - /Aligned.sortedByCoord.out/
    """
    picard MarkDuplicates \\
        INPUT=$bam \\
        OUTPUT=${sample}.markDups.bam \\
        METRICS_FILE=${sample}.markDups_metrics.txt \\
        REMOVE_DUPLICATES=false \\
        ASSUME_SORTED=true \\
        PROGRAM_RECORD_ID='null' \\
        VALIDATION_STRINGENCY=LENIENT \\
	CREATE_INDEX=TRUE
    """
}


/*
 * Add Read Groups
 */
process addRG {
    tag "${sample}"
    publishDir "${params.outdir}/markDuplicates", mode: 'copy'

    input:
    set val(sample), file(markdup_bam), file(markdup_bam_ind) from bam_md

    output:
    set val("$sample"), file("${sample}.markDups.wRG.bam"), file("${sample}.markDups.wRG.bai") into bam_md_rg

    script:
    """
    picard AddOrReplaceReadGroups \\
        INPUT=$markdup_bam \\
        OUTPUT=${sample}.markDups.wRG.bam \\
        RGLB=LB RGPL=ILLUMINA RGPU=PU RGSM=${sample} \\
	CREATE_INDEX=TRUE
    """
}



/*
 * Split and trim
 */
process splitandtrim {
    tag "${sample}"
    publishDir "${params.outdir}/splitandtrim", mode: 'copy'

    input:
    set val(sample), file(markdup_bam), file(markdup_bam_ind) from bam_md_rg

    output:
    set val("$sample"), file("${sample}_split.bam"), file("${sample}_split.bai") into bam_split

    script:
    """
    gatk SplitNCigarReads \\
	--TMP_DIR tmp \\
	-I $markdup_bam \\
	-R $params.fasta \\
	-O ${sample}_split.bam \\
    	-OBI 
    """
}


/*
 * Recalibrate BAM file with known variants and BaseRecalibrator
 *
*/

process recalibrate {
    tag "${sample}"
    publishDir "${params.outdir}/recalibration", mode: 'copy'
	
    input:
    set val(sample), file(split_bam), file(split_bam_ind) from bam_split

    output:
    set val(sample), file("${sample}_recal.bam"), file("${sample}_recal.bai") into bam_vcall
    file '.command.log' into gatk_base_recalibration_results
	
    script:
    """
    gatk BaseRecalibrator \\
        -I $split_bam \\
        -R $params.fasta \\
        -O ${sample}_table.recal \\
	--known-sites $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g

    gatk ApplyBQSR \\
        -R $params.fasta \\
        -I $split_bam \\
        --bqsr-recal-file ${sample}_table.recal \\
        -O ${sample}_recal.bam \\
        --create-output-bam-index true \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}


// Call variants
process variantCall {
    tag "$sample"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    set val(sample), file(recal_bam), file(recal_bam_ind) from bam_vcall

    output:
    set val(sample), file("${sample}_raw_variants.vcf"), file("${sample}_raw_variants.vcf.idx") into raw_variants_gvcf

    script:
    """
    gatk HaplotypeCaller \\
        -I $recal_bam \\
        -R $params.fasta \\
        -O ${sample}_raw_variants.vcf \\
	--dont-use-soft-clipped-bases \\
	-stand-call-conf 20.0 \\
        -ERC GVCF \\
        -OBI \\
	--annotation MappingQualityRankSumTest \\
        --annotation QualByDepth \\
        --annotation ReadPosRankSumTest \\
        --annotation RMSMappingQuality \\
        --annotation FisherStrand \\
        --annotation Coverage \\
        --dbsnp $params.dbsnp \\
        --verbosity INFO \\
        --java-options -Xmx${task.memory.toGiga()}g
    """
}

/*
 * Genotype generate GVCFs using GATK's GenotypeGVCFs
 * 
*/ 

process genotypegvcfs{
    tag "$sample"
    publishDir "${params.outdir}/variants", mode: 'copy' 

    input:
    set val(sample), file(raw_vcf), file(raw_vcf_idx) from raw_variants_gvcf

    output:
    set val(sample), file("${sample}_gvcf.vcf"), file("${sample}_gvcf.vcf.idx") into raw_gvcfs

    script:
    """
    gatk GenotypeGVCFs \\
    -R $params.fasta \\
    --dbsnp $params.dbsnp \\
    -V $raw_vcf \\
    -O ${sample}_gvcf.vcf 
    """
}


// Hard-filter variants
process filterVariants {
    tag "$sample"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    set val(sample), file(raw_gvcf), file(raw_gvcf_idx) from raw_gvcfs

    output:
    set val(sample), file("${sample}_filtered.vcf"), file("${sample}_filtered.vcf.idx") into filtered_variants

    script:
    """
    gatk VariantFiltration \\
        -R $params.fasta \\
        -V $raw_gvcf \\
        -O ${sample}_filtered.vcf \\
        -window 35 -cluster 3 \\
	--filter-name FS \\
	--filter-expression "FS > 30.0" \\
	--filter-name QD \\
	--filter-expression "QD < 2.0"
    """
}


// Annotate with annovar
process annotateVariants {
    tag "$sample"
    publishDir "${params.outdir}/variants", mode: 'copy'

    input:
    set val(sample), file(filtered_vcf), file(filtered_vcf_idx) from filtered_variants

    output:
    set val(sample), file("${sample}.hg38_multianno.vcf"), file("${sample}.hg38_multianno.txt") into annotated_variants

    script:
    """
    table_annovar.pl ${filtered_vcf} $params.annovar_anno \\
    --outfile ${sample} \\
    --buildver $params.genome \\
    --protocol refGene,cosmic87_coding,cosmic87_noncoding,clinvar_20180603,avsnp150,1000g2015aug_all,gnomad_genome,dbnsfp35a,dbscsnv11 \\
    --operation g,f,f,f,f,f,f,f,f \\
    --vcfinput \\
    --polish \\
    --remove
    """
}



/*
 * Parse software version numbers
 */
process get_software_versions {

    output:
    file 'software_versions_mqc.yaml' into software_versions_yaml

    script:
    """
    echo $workflow.manifest.version > v_pipeline.txt
    echo $workflow.nextflow.version &> v_nextflow.txt
    fastqc --version &> v_fastqc.txt
    cutadapt --version &> v_cutadapt.txt
    trim_galore --version &> v_trim_galore.txt
    STAR --version &> v_star.txt
    picard MarkDuplicates --version &> v_markduplicates.txt  || true
    gatk ApplyBQSR --help 2>&1 | grep Version: > v_gatk.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt
    scrape_software_versions.py &> software_versions_mqc.yaml
    """
}


/*
 * Pipeline parameters to go into MultiQC report
 */
process workflow_summary_mqc {

    output:
    file 'workflow_summary_mqc.yaml' into workflow_summary_yaml

    exec:
    def yaml_file = task.workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'ucl-blic-rnaseq-variant-calling-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'UCL-BLIC/rnaseq_variant_calling Workflow Summary'
    section_href: 'https://github.com/UCL-BLIC/rnaseq_variant_calling'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k,v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()
}


/*
 * STEP 12 MultiQC
 */
process multiqc {
    publishDir "${params.outdir}/MultiQC", mode: 'copy'

    input:
    file multiqc_config
    file ('fastqc/*') from fastqc_results.collect().ifEmpty([])
    file ('trimgalore/*') from trimgalore_results.collect()
    file ('alignment/*') from alignment_logs.collect()
    file ('gatk_base_recalibration/T*') from gatk_base_recalibration_results.collect()
    file ('software_versions/*') from software_versions_yaml.collect()
    file ('workflow_summary/*') from workflow_summary_yaml.collect()

    output:
    file "*multiqc_report.html" into multiqc_report
    file "*_data"

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    """
    multiqc . -f $rtitle $rfilename --config $multiqc_config \\
        -m custom_content -m picard -m star -m cutadapt -m fastqc
    """
}

/*
 * STEP 13 - Output Description HTML
 */
process output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: 'copy'

    input:
    file output_docs

    output:
    file "results_description.html"

    script:
    """
    markdown_to_html.r $output_docs results_description.html
    """
}


/*
 * Completion e-mail notification
 */
workflow.onComplete {

//////////////////////////////////////////////////////////////////////////////////////
    if(skipped_poor_alignment.size() > 0){
        log.info "[nfcore/rnaseq] WARNING - ${skipped_poor_alignment.size()} samples skipped due to poor alignment scores!"
    }
//////////////////////////////////////////////////////////////////////////////////////

    log.info "[nfcore/rnaseq] Pipeline Complete"

}

