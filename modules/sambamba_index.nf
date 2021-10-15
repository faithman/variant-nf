/*

Sambamba index module

*/

process SAMBAMBA_INDEX {
	
	tag { sm }

	publishDir "${params.output}/bams", mode: "copy"

	label "large"

	input:
	tuple val(sm), path(bam)

	output:
	tuple val(sm), path(bam), path("*.bai"), emit: bam_indexed

	"""
	sambamba \\
		index ${bam}
	"""
}