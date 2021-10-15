/*

Multiqc

*/

process MULTIQC_CALLING {

  tag { "Multiqc plot" }

  label 'large'

  publishDir "${params.output}/report", mode: 'copy'

  input:
      path("soft-filter.stats.txt")
      path("hard-filter.stats.txt")

  output:
      file("multiqc.html")
      file("multiqc_data/*.json")

  """
      multiqc -k json --filename multiqc.html .
  """

}
