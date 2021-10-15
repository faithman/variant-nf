/*

Picard markduplicates module

*/

process PICARD_MARKDUPLICATES {
	
	tag { sm }

	//publishDir "${params.output}/bams", mode: "copy", pattern: "*.bam"
	publishDir "${params.output}/metrics", mode: "copy", pattern: "*.metrics.txt"

	label "large"

	input:
	//tuple val(row.sm), path(bam), path(bai)
	tuple val(sm), val(replicate), path(bam), path(bai)

	output:
	tuple val(sm), path("*.mkdup.bam"), emit: bam_marked
	path("*.metrics.txt"), emit: metrics

	"""
	picard \\
		-Xmx64g \\
		-Djava.io.tmpdir=/home/wangye/tmp \\
		MarkDuplicates \\
		INPUT=${bam} \\
		OUTPUT=${sm}_${replicate}.mkdup.bam \\
		METRICS_FILE=${sm}.MarkDuplicates.metrics.txt
	"""
}
