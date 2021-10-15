/*

Multiqc

*/

process MULTIQC_MAPPING {

  tag { "Multiqc plot" }

  label 'large'

  publishDir "${params.output}/report", mode: 'copy'

  input:
      path("*_fastp.json")
      path("*.metrics.txt")
      path("*stats")
      path("*idxstats")
      path("*flagstat")
      path("*.mosdepth.summary.txt")

  output:
      file("multiqc.html")
      file("multiqc_data/*.json")

  """
      multiqc -k json --filename multiqc.html .
  """

}
