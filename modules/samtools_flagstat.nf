/* 

Samtools flagstat module

*/

process SAMTOOLS_FLAGSTAT {

	tag {sm}

	label 'large'

	publishDir "${params.output}/stats", mode: "copy", pattern: "*.flagstat"

	input:
		tuple val(sm), path(bam), path(bai)

	output:
		path("${bam}.flagstat"), emit: flagstat

	"""
		samtools flagstat ${bam} > ${bam}.flagstat
	"""
}