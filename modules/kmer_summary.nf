/*

Kmer summary module

*/

process KMER_SUMMARY {

    tag "kmer_summary"

    label 'small'

    publishDir "${params.output}/kmer", mode: 'copy'

    input:
        file("*kmer.tsv")

    output:
        file("kmers.tsv")

    """
        cat <(echo "kmer\tfrequency\tgroup\tsm\treads ") *.tsv > kmers.tsv
    """

}
