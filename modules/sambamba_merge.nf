/*

Sambamba merge bam module

*/


process SAMBAMBA_MERGE {

    tag { sm }

    label 'large'

    publishDir "${params.output}/bams/merged", mode: "copy"

    input:
        tuple val(sm), path(bam), path(bai), val(n_count)

    output:
        tuple val(sm), path("${sm}.bam"), path("${sm}.bam.bai"), emit: merged_bam

    script:
        if (n_count == 1)

            """
                mv ${bam} ${sm}.bam
                mv ${bai} ${sm}.bam.bai
            """
        else

            """
                sambamba merge --nthreads=${task.cpus} --show-progress ${sm}.bam ${bam}
                sambamba index --nthreads=${task.cpus} ${sm}.bam
            """
}
