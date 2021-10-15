#!/usr/bin/env nextflow

/*
    variant-nf Pipeline
    Chengdu Research Base of Giant Panda Breeding
    Authors:
    - Ye Wang <yewangfaith@gmail.com>
*/
nextflow.enable.dsl=2

// This pipeline requires NXF_VER 20.01.0-rc1 or later

/*
    Params
*/


date = new Date().format( 'yyyyMMdd' )
params.time = "${date}"
params.debug =  false
params.email = "yewang_faith@hotmail.com"

// Check that reference exists
params.genome_dir = ""
params.reference_genome = ""
params.gtf = ""
params.cpus = ""
params.memory = ""
params.mapping_only = ""
params.calling_only = ""

if (params.calling_only.toString() == "true" && params.debug.toString() == "") {

    params.bam_sheet = "${workflow.projectDir}/bam_list.tsv"

}


//reference = file(params.reference, checkIfExists: true)

// Debug
if (params.debug.toString() == "true") {
    params.output = "debug-${date}"
    params.sample_sheet = "${workflow.projectDir}/test_data/test_sample_sheet.tsv"
    params.bam_sheet = "test_data/test_bam_list.tsv"
} else {
    // The strain sheet that used for 'production' is located in the root of the git repo
    params.output = "standard-${date}"
    params.sample_sheet = "${workflow.projectDir}/sample_sheet.tsv"
}


def log_summary() {

    out =  '''

▗  ▖ ▗▖ ▗▄▄ ▗▄▄  ▗▖ ▗▖ ▖▄▄▄▖    ▗▖ ▖▗▄▄▖
▝▖▗▘ ▐▌ ▐ ▝▌ ▐   ▐▌ ▐▚ ▌ ▐      ▐▚ ▌▐
 ▌▐  ▌▐ ▐▄▄▘ ▐   ▌▐ ▐▐▖▌ ▐      ▐▐▖▌▐▄▄▖
 ▚▞  ▙▟ ▐ ▝▖ ▐   ▙▟ ▐ ▌▌ ▐   ▀▘ ▐ ▌▌▐
 ▐▌ ▐  ▌▐  ▘▗▟▄ ▐  ▌▐ ▐▌ ▐      ▐ ▐▌▐
                                              
'''

out += """
    parameters               description                 Set/Default
    ==========               ===========                 ========================
    Debug                    Debug or not                ${params.debug}
    output                   Release Directory           ${params.output}
    bam_list		     list for bam file		 ${params.bam_sheet}
    sample_sheet             sample sheet                ${params.sample_sheet}
    genome_dir               Genome and index            ${params.genome_dir}
    username                                             ${"whoami".execute().in.text}

    Nextflow Run
    ---------------
    ${workflow.commandLine}
    run name                                             ${workflow.runName}
    scriptID                                             ${workflow.scriptId}
    git commit                                           ${workflow.commitId}
    container                                            ${workflow.container}
---
"""
out
}

log.info(log_summary())


// Include processes from modules
include { FASTP_TRIM } from './modules/fastp_trim.nf' params(params)
include { BWA_MEM } from './modules/bwa.nf' params(params)
include { PICARD_MARKDUPLICATES } from './modules/picard_markdup.nf' params(params)
include { SAMBAMBA_INDEX } from './modules/sambamba_index.nf' params(params)
//include { SAMBAMBA_INDEX } as { SAMBAMBA_INDEX_MERGED } from './modules/sambamba_index.nf' params(params)
include { SAMTOOLS_STATS } from './modules/samtools_stats.nf' params(params)
include { SAMTOOLS_IDXSTATS } from './modules/samtools_idxstats.nf' params(params)
include { SAMTOOLS_FLAGSTAT } from './modules/samtools_flagstat.nf' params(params)
include { MOSDEPTH } from './modules/mosdepth.nf' params(params)
include { KMER_COUNTING } from './modules/kmer_counting.nf' params(params)
include { KMER_SUMMARY } from './modules/kmer_summary.nf' params(params)
include { MOSDEPTH_PLOT } from './modules/mosdepth_plot.nf' params(params)
include { SAMBAMBA_MERGE } from './modules/sambamba_merge.nf' params(params)
include { GET_CONTIGS } from './modules/gatk_short_variant_calling.nf' params(params)
include { HAPLOTYPECALLER_INDIVIDUAL } from './modules/gatk_short_variant_calling.nf' params(params)
include { CONCAT_INDIVIDUAL_GVCFS } from './modules/gatk_short_variant_calling.nf' params(params)
include { IMPORT_GENOMICS_DB } from './modules/gatk_short_variant_calling.nf' params(params)
include { GENOTYPE_COHORT_GVCF } from './modules/gatk_short_variant_calling.nf' params(params)
include { CONCATENATE_VCF } from './modules/gatk_short_variant_calling.nf' params(params)
include { SOFT_FILTER } from './modules/gatk_short_variant_calling.nf' params(params)
include { HARD_FILTER } from './modules/gatk_short_variant_calling.nf' params(params)
include { MULTIQC } from './modules/multiqc.nf' params(params)
include { MULTIQC_MAPPING } from './modules/multiqc_mapping.nf' params(params)
include { MULTIQC_CALLING } from './modules/multiqc_calling.nf' params(params)


if (workflow.profile == "") {
    println "Must set -profile: debug, standard"
    exit 1
}

// Read sample sheet into channel
sample_sheet = Channel 
    .fromPath(params.sample_sheet, checkIfExists: true)
    .ifEmpty {exit 1, "sample sheet not found"}
    .splitCsv(header:true, sep: "\t")
    .map { row -> row.fq1 = row.fq1;row}
    .map { row -> row.fq2 = row.fq2;row}
    .map { row -> [row, file(row.fq1), file(row.fq2)]}
    //.view()

if (params.calling_only.toString() == "true") {

    bam_list = Channel 
    .fromPath(params.bam_sheet, checkIfExists: true)
    .ifEmpty {exit 1, "sample sheet not found"}
    .splitCsv(header:true, sep: "\t")
    .map { row -> row.bam = row.bam;row}
    .map { row -> row.bai = row.bai;row}
    .map { row -> [row.sm, file(row.bam), file(row.bai)]}

}


//Read reference genome into channel

//reference_genome = Channel.fromPath(params.genome_dir)
    
// Workflow

workflow MAPPING {

    take: sample

    main:

    // Trim raw fastq files
    FASTP_TRIM( sample )

    // Perform mapping
    BWA_MEM( FASTP_TRIM.out.fastq_trimmed )

    // K-mer counting
    KMER_COUNTING( FASTP_TRIM.out.fastq_trimmed )

    // K-mer summary
    KMER_SUMMARY( KMER_COUNTING.out.kmer_counting)

    // Mark dup
    PICARD_MARKDUPLICATES( BWA_MEM.out.bwa )

    // Index marked BAM
    SAMBAMBA_INDEX( PICARD_MARKDUPLICATES.out.bam_marked )

    // Merge bams from one biological sample
    merge_bam = SAMBAMBA_INDEX.out.bam_indexed.map { sm, bam, bai -> [sm, bam, bai] }
                        .groupTuple( by:0 )
                        .map { sm, bam, bai -> [sm, bam, bai, bam.size()] }

    SAMBAMBA_MERGE( merge_bam )

    // BAM stats
    SAMTOOLS_STATS( SAMBAMBA_MERGE.out.merged_bam )

    // BAM idxstats
    SAMTOOLS_IDXSTATS( SAMBAMBA_MERGE.out.merged_bam )

    // Flagstats
    SAMTOOLS_FLAGSTAT( SAMBAMBA_MERGE.out.merged_bam )

    // BAM coverage
    MOSDEPTH( SAMBAMBA_MERGE.out.merged_bam )

    // Plot coverage
    MOSDEPTH_PLOT( MOSDEPTH.out.mosdepth_dist.collect() )

    emit:

    bams = SAMBAMBA_MERGE.out.merged_bam
    fastp = FASTP_TRIM.out.fastp_json.collect()
    picard_mkdup = PICARD_MARKDUPLICATES.out.metrics.collect()
    samtools_stats = SAMTOOLS_STATS.out.stats
    samtools_idx = SAMTOOLS_IDXSTATS.out.idxstats
    samtools_flag = SAMTOOLS_FLAGSTAT.out.flagstat
    mosdepth = MOSDEPTH.out.mosdepth_summary.collect()

}


workflow VARIANT_CALLING {

    take: bam_file

    main:

    // Get contigs from bam file
    GET_CONTIGS( bam_file.first() )

    // Split contigs
    contigs = GET_CONTIGS.out.contigs.splitText() { it.strip() }

    // Merge bam list with contigs
    new_sample_sheet = bam_file.groupTuple().combine( contigs )
    
    // Call variants for each sample and each contig
    HAPLOTYPECALLER_INDIVIDUAL( new_sample_sheet )

    // Generate gvcf list
    gvcf_list = HAPLOTYPECALLER_INDIVIDUAL.out.individual_contig_gvcfs
                                              .groupTuple()
                                              .map { sm, gvcf -> [sm, gvcf]}
                                              .combine( GET_CONTIGS.out.contigs)
                                              

    // Concat gvcf files for each sample
    CONCAT_INDIVIDUAL_GVCFS( gvcf_list )

    // generate sample map
    sample_map = bam_file.groupTuple().map { "${it[0]}\t${it[0]}.g.vcf.gz" }.collectFile(name: "sample_map.tsv", newLine: true)

    // db_list
    dblist = CONCAT_INDIVIDUAL_GVCFS.out.ind_gvcfs
                                        .flatten()
                                        .toList()
                                        .map { [it] }
                                        .combine( contigs )
                                        .combine( sample_map )
                                        
    // Import genomics db
    IMPORT_GENOMICS_DB( dblist )

    // Genotype gVCFs
    GENOTYPE_COHORT_GVCF( IMPORT_GENOMICS_DB.out.ind_db )

    // Concatenate all vcf files
    CONCATENATE_VCF( GENOTYPE_COHORT_GVCF.out.cohort_vcf.collect() )

    // Perform soft filter
    SOFT_FILTER( CONCATENATE_VCF.out.concatenated_vcf )

    // Perform Soft filter
    hard_filter_input = SOFT_FILTER.out.soft_filter_vcf
                                       .combine( GET_CONTIGS.out.contigs )

    HARD_FILTER( hard_filter_input )

    emit:

    soft_filter = SOFT_FILTER.out.soft_vcf_stats
    hard_filter = HARD_FILTER.out.hard_vcf_stats

}


workflow {

    if ( params.mapping_only.toString() == "true" && params.calling_only == "" ) {
        MAPPING( sample_sheet )

        MULTIQC_MAPPING (
            MAPPING.out.fastp, 
            MAPPING.out.picard_mkdup,
            MAPPING.out.samtools_stats,
            MAPPING.out.samtools_idx,
            MAPPING.out.samtools_flag,
            MAPPING.out.mosdepth )

    } else if ( params.mapping_only == "" && params.calling_only.toString() == "true" ) {

        VARIANT_CALLING( bam_list )

        MULTIQC_CALLING (
            VARIANT_CALLING.out.soft_filter, 
            VARIANT_CALLING.out.hard_filter )

    } else {

        MAPPING( sample_sheet )

        VARIANT_CALLING( MAPPING.out.bams )

        MULTIQC (
            MAPPING.out.fastp, 
            MAPPING.out.picard_mkdup,
            MAPPING.out.samtools_stats,
            MAPPING.out.samtools_idx,
            MAPPING.out.samtools_flag,
            MAPPING.out.mosdepth,
            VARIANT_CALLING.out.soft_filter, 
            VARIANT_CALLING.out.hard_filter )
    }
}
    /*

*/

/*
workflow {

    // Generate a summary of the current run
    summary(Channel.from("run"))

    /*=====================
      Short reads mapping
    =======================*/
/*
    // Trim raw fastq files
    FASTP_TRIM( sample_sheet )

    // Perform mapping
    BWA_MEM( FASTP_TRIM.out.fastq_trimmed )

    // K-mer counting
    KMER_COUNTING( FASTP_TRIM.out.fastq_trimmed )

    // K-mer summary
    KMER_SUMMARY( KMER_COUNTING.out.kmer_counting)

    // Mark dup
    PICARD_MARKDUPLICATES( BWA_MEM.out.bwa )

    // Index marked BAM
    SAMBAMBA_INDEX( PICARD_MARKDUPLICATES.out.bam_marked )

    // Merge bams from one biological sample
    merge_bam = SAMBAMBA_INDEX.out.bam_indexed.map { sm, bam, bai -> [sm, bam, bai] }
                        .groupTuple( by:0 )
                        .map { sm, bam, bai -> [sm, bam, bai, bam.size()] }

    SAMBAMBA_MERGE( merge_bam )

    // BAM stats
    SAMTOOLS_STATS( SAMBAMBA_MERGE.out.merged_bam )

    // BAM idxstats
    SAMTOOLS_IDXSTATS( SAMBAMBA_MERGE.out.merged_bam )

    // Flagstats
    SAMTOOLS_FLAGSTAT( SAMBAMBA_MERGE.out.merged_bam )

    // BAM coverage
    MOSDEPTH( SAMBAMBA_MERGE.out.merged_bam )

    // Plot coverage
    MOSDEPTH_PLOT( MOSDEPTH.out.mosdepth_dist.collect() )

    /*=====================
      Short variant calling
    =======================*/
/*
    // Get contigs from bam file
    GET_CONTIGS( SAMBAMBA_MERGE.out.merged_bam.first() )

    // Split contigs
    contigs = GET_CONTIGS.out.contigs.splitText() { it.strip() }

    // Merge bam list with contigs
    new_sample_sheet = SAMBAMBA_MERGE.out.merged_bam.groupTuple().combine( contigs )
    
    // Call variants for each sample and each contig
    HAPLOTYPECALLER_INDIVIDUAL( new_sample_sheet )

    // Generate gvcf list
    gvcf_list = HAPLOTYPECALLER_INDIVIDUAL.out.individual_contig_gvcfs
                                              .groupTuple()
                                              .map { sm, gvcf -> [sm, gvcf]}
                                              .combine( GET_CONTIGS.out.contigs)
                                              

    // Concat gvcf files for each sample
    CONCAT_INDIVIDUAL_GVCFS( gvcf_list )

    // generate sample map
    sample_map = SAMBAMBA_MERGE.out.merged_bam.groupTuple().map { "${it[0]}\t${it[0]}.g.vcf.gz" }.collectFile(name: "sample_map.tsv", newLine: true)

    // db_list
    dblist = CONCAT_INDIVIDUAL_GVCFS.out.ind_gvcfs
                                        .flatten()
                                        .toList()
                                        .map { [it] }
                                        .combine( contigs )
                                        .combine( sample_map )
                                        
    // Import genomics db
    IMPORT_GENOMICS_DB( dblist )

    // Genotype gVCFs
    GENOTYPE_COHORT_GVCF( IMPORT_GENOMICS_DB.out.ind_db )

    // Concatenate all vcf files
    CONCATENATE_VCF( GENOTYPE_COHORT_GVCF.out.cohort_vcf.collect() )

    // Perform soft filter
    SOFT_FILTER( CONCATENATE_VCF.out.concatenated_vcf )

    // Perform Soft filter
    hard_filter_input = SOFT_FILTER.out.soft_filter_vcf
                                       .combine( GET_CONTIGS.out.contigs )

    HARD_FILTER( hard_filter_input )

    // Generate multqic report
    MULTIQC (FASTP_TRIM.out.fastp_json.collect(), 
        PICARD_MARKDUPLICATES.out.metrics.collect(),
        SAMTOOLS_STATS.out.stats,
        SAMTOOLS_IDXSTATS.out.idxstats,
        SAMTOOLS_FLAGSTAT.out.flagstat,
        MOSDEPTH.out.mosdepth_summary.collect(),
        SOFT_FILTER.out.soft_vcf_stats,
        HARD_FILTER.out.hard_vcf_stats
        )
}
*/

process summary {

    // This process is for a summary of current run.
    
    executor 'local'

    publishDir "${params.output}", mode: 'copy'
    
    input:
        val(run)

    output:
        path("sample_sheet.tsv")
        path("summary.txt")

    """
        echo '''${log_summary()}''' > summary.txt
        cat ${params.sample_sheet} > sample_sheet.tsv
    """

}

workflow.onComplete {

    def msg = """\
        Pipeline execution summary
        ---------------------------
        Completed at: ${workflow.complete}
        Duration    : ${workflow.duration}
        Success     : ${workflow.success}
        workDir     : ${workflow.workDir}
        exit status : ${workflow.exitStatus}
        """
        .stripIndent()

    sendMail(to: params.email, subject: 'My pipeline execution', body: msg)
}
