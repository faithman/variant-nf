/*

Samtool stats module

*/

process SAMTOOLS_STATS {

    tag {sm}

    publishDir "${params.output}/stats", mode: "copy", pattern: "*.stats"

    label "large"

    input:
    tuple val(sm), path(bam), path(bai)
    
    output:
    path("*.stats"), emit: stats


    """
    samtools stats ${bam} > ${bam}.stats
    """
}