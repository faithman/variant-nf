/*
 
    samtools idx stats

*/

process SAMTOOLS_IDXSTATS {
    
    tag {sm}

    label 'large'

    publishDir "${params.output}/stats", mode: "copy", pattern: "*.idxstats"

    input:
        tuple val(sm), path(bam), path(bai)

    output:
        path("*.idxstats"), emit: idxstats

    """
        samtools idxstats ${bam} > ${bam}.idxstats
    """
}
