/*

Coverage

*/

process MOSDEPTH {
	
	tag {sm}

	label 'large'

	publishDir "${params.output}/coverage", mode: 'copy', pattern: "${sm}.*"

	input:
		tuple val(sm), path(bam), path(bai)

	output:
		path("${sm}.mosdepth.summary.txt"), emit: mosdepth_summary
		path("${sm}.mosdepth.global.dist.txt"), emit: mosdepth_dist
		path("${sm}.per-base.bed.gz"), emit: mosdepth_bed
		path("${sm}.per-base.bed.gz.csi"), emit: mosdepth_csi

	"""
		export MOSDEPTH_PRECISION=5
		mosdepth --threads ${task.cpus} \\
				 ${sm} \\
				 ${bam}
	"""
}