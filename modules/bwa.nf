/*

BWA-MEM mapping module

*/

process BWA_MEM {

    tag { row.sm }
    
    label 'large'

    input:
        tuple val(row), path(fq1), path(fq2)
        //path reference_genome


    output:
        tuple val(row.sm), val(row.replicate), path("${row.sm}.bam"), path("${row.sm}.bam.bai"), emit: bwa
	    tuple val(row.sm), path("${row.sm}.bam"), emit: bwa_nonindex

	script:
		// Construct read group
		RG = ["@RG",
			  "ID:${row.sm}",
			  "SM:${row.group}"].join("\\t")

    """
        bwa mem -t ${task.cpus} -R '${RG}' ${params.genome_dir} ${fq1} ${fq2} | \\
        sambamba view --nthreads=${task.cpus} --show-progress --sam-input --format=bam --with-header /dev/stdin | \\
        sambamba sort --nthreads=${task.cpus} --show-progress --tmpdir=. --out=${row.sm}.bam /dev/stdin
        sambamba index --nthreads=${task.cpus} ${row.sm}.bam
        
        if [[ ! \$(samtools view ${row.sm}.bam | head -n 10) ]]; then
            exit 1;
        fi
    """
}

