/*

picard validatesamfile module

*/

process validatebam {

    tag {sm}

    label 'large'

    publishDir "${params.output}/bam_validate", mode: "copy", pattern: "*.validatesamfile.txt"

    input:
        tuple val(sm), path(bam), file(bai)

    output:
        tuple sm, path("*.validatesamfile.txt"), emit: bam_validate

    """
        picard ValidateSamFile I=${bam} MODE=SUMMARY > ${bam}.validatesamfile.txt
    """
}
