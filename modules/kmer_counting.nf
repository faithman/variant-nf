/*

Fastq-tools kmer module

*/

process KMER_COUNTING {

    tag {row.sm}
 
    label "large"

    input:
        tuple val(row), path(fq1), path(fq2)

    output:
        tuple val(row.sm), file("${row.sm}.kmer.tsv"), emit: kmer_counting

    """
        export OFS="\t"

        fq_wc=`zcat ${fq1} | awk 'NR % 4 == 0' | wc -l`

        zcat ${fq1} ${fq2} | \\
        fastq-kmers -k 6 | \\
        awk -v OFS="\t" -v ID=${row.sm} -v SM=${row.group} -v fq_wc="\${fq_wc}" 'NR > 1 { print \$0, SM, ID, fq_wc }' - > ${row.sm}.kmer.tsv
    """
}
