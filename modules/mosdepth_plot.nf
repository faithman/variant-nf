/*

Plot coverage

*/

process MOSDEPTH_PLOT {

  tag { "mosdepth_plot" }

  publishDir "${params.output}/coverage", mode: 'copy'

  input:
    file("*.dist.txt") 

  output:
    file("dist.html")

  script:
  """
    plot-dist.py *.txt 
  """

}