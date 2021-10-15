/*

Get contigs 

*/

process GET_CONTIGS {

  tag { sm }

  label 'small'

  input:
      tuple val(sm), path(bam), path(bai)

  output:
      path("contigs.txt"), emit: contigs

  """
      samtools idxstats ${bam} | cut -f 1 | grep -v "*" > contigs.txt
  """
}

/*

Run gatk haplotypecaller to call each file at contig level

*/

process HAPLOTYPECALLER_INDIVIDUAL {

  tag { "${sm}:${contig_number}" }

  label 'large'

  input:
      tuple val(sm), path(bam), path(bai), val(contig_number)

  output:
      tuple val(sm), path("${contig_number}.g.vcf.gz"), emit: individual_contig_gvcfs

  """
      gatk HaplotypeCaller --java-options "-Xmx32g -Xms8g -XX:ConcGCThreads=${task.cpus}" \\
        --emit-ref-confidence GVCF \\
        --annotation DepthPerAlleleBySample \\
        --annotation Coverage \\
        --annotation GenotypeSummaries \\
        --annotation TandemRepeat \\
        --annotation StrandBiasBySample \\
        --annotation ChromosomeCounts \\
        --annotation ReadPosRankSumTest \\
        --annotation AS_ReadPosRankSumTest \\
        --annotation AS_QualByDepth \\
        --annotation AS_StrandOddsRatio \\
        --annotation AS_MappingQualityRankSumTest \\
        --annotation DepthPerSampleHC \\
        --annotation-group StandardAnnotation \\
        --annotation-group AS_StandardAnnotation \\
        --annotation-group StandardHCAnnotation \\
        --do-not-run-physical-phasing \\
        -R ${params.reference_genome} \\
        -I ${bam} \\
        -L ${contig_number} \\
        -O ${contig_number}.g.vcf

      bcftools view -O z ${contig_number}.g.vcf > ${contig_number}.g.vcf.gz

      rm ${contig_number}.g.vcf
  """
}

/*

Merge individual gVCF files

*/

process CONCAT_INDIVIDUAL_GVCFS {

  tag { sm }

  label 'large'

  input:
      tuple val(sm), path("*"), path(contigs)

  output:
      tuple path("${sm}.g.vcf.gz"), path("${sm}.g.vcf.gz.tbi"), emit: ind_gvcfs

  """
      awk '{ print \$0 ".g.vcf.gz" }' ${contigs} > contig_set.tsv
      bcftools concat -O z --file-list contig_set.tsv > ${sm}.g.vcf.gz
      bcftools index --tbi ${sm}.g.vcf.gz
  """
}

/*

Import genomics database

*/

  process IMPORT_GENOMICS_DB {

    tag { contig }
    
    label 'large'

    input:
        tuple path(vcfs), val(contig), path(sample_map)

    output:
        tuple val(contig), file("${contig}.db"), emit: ind_db

    """
        gatk --java-options "-Xmx64g -Xms8g -XX:ConcGCThreads=${task.cpus}" \\
             GenomicsDBImport --genomicsdb-workspace-path ${contig}.db \\
                              --batch-size 12 \\
                              -L ${contig} \\
                              --sample-name-map ${sample_map} \\
                              --reader-threads ${task.cpus}
    """

}

/*

Genotype gVCFs

*/

process GENOTYPE_COHORT_GVCF {

  tag { contig }

  label 'large'

  input:
      tuple val(contig), file("${contig}.db")

  output:
      tuple file("${contig}_cohort.vcf.gz"), file("${contig}_cohort.vcf.gz.csi"), emit: cohort_vcf

  """
      gatk --java-options "-Xmx64g -Xms8g -XX:ConcGCThreads=${task.cpus}" \\
           GenotypeGVCFs -R ${params.reference_genome} \\
                         -V gendb://${contig}.db \\
                         -G StandardAnnotation \\
                         -G AS_StandardAnnotation \\
                         -G StandardHCAnnotation \\
                         -L ${contig} \\
                         -O ${contig}_cohort.vcf

      bcftools view -O z --min-af 0.00001 --threads=${task.cpus} ${contig}_cohort.vcf > ${contig}_cohort.vcf.gz
      bcftools index ${contig}_cohort.vcf.gz
  """
}

/*

Genotype gVCFs

*/

process CONCATENATE_VCF {

  tag { "concatenate_vcf" }

  label 'large'

  input:
      path("*")

  output:
      tuple path("raw.vcf.gz"), path("raw.vcf.gz.tbi"), emit: concatenated_vcf

  """
      ls *_cohort.vcf.gz > contig_set.tsv
      bcftools concat -O z --file-list contig_set.tsv > raw.vcf.gz 
      bcftools index --tbi raw.vcf.gz
  """
}

/*

Soft filter

*/

process SOFT_FILTER {

  tag { "Performing soft filter" }

  label 'large'

  publishDir "${params.output}/variation/soft_filter", mode: 'copy'

  input:
      tuple path("raw.vcf.gz"), path("raw.vcf.gz.tbi")

  output:
      tuple path("${params.time}.soft-filter.vcf.gz"), path("${params.time}.soft-filter.vcf.gz.tbi"), emit: soft_filter_vcf
      path("${params.time}.soft-filter.vcf.gz.csi")
      path("${params.time}.soft-filter.stats.txt"), emit: soft_vcf_stats
      path("filter_stats.tsv")

  """
      gatk --java-options "-Xmx64g -Xms8g -XX:ConcGCThreads=${task.cpus}" \\
          VariantFiltration \\
              -R ${params.reference_genome} \\
              --variant raw.vcf.gz \\
              --genotype-filter-expression "DP < ${params.min_depth}" --genotype-filter-name "DP_min_depth" \\
              --filter-expression "QUAL < ${params.qual}" --filter-name "QUAL_quality" \\
              --filter-expression "FS > ${params.fisherstrand}" --filter-name "FS_fisherstrand" \\
              --filter-expression "QD < ${params.quality_by_depth}" --filter-name "QD_quality_by_depth" \\
              --filter-expression "SOR > ${params.strand_odds_ratio}" --filter-name "SOR_strand_odds_ratio" \\
              -O out.vcf

      bcftools view -O z out.vcf > out.vcf.gz
      bcftools index --tbi out.vcf.gz

      # filtering high missing calls
      bcftools filter --threads ${task.cpus} --soft-filter='high_missing' --mode + --include 'F_MISSING <= ${params.high_missing}' out.vcf.gz -O z > ${params.time}.soft-filter.vcf.gz

      bcftools index --tbi ${params.time}.soft-filter.vcf.gz
      bcftools index ${params.time}.soft-filter.vcf.gz

      rm out.vcf out.vcf.gz out.vcf.gz.tbi

      bcftools stats --threads ${task.cpus} -s- --verbose ${params.time}.soft-filter.vcf.gz > ${params.time}.soft-filter.stats.txt

      {
        echo -e 'QUAL\\tQD\\tSOR\\tFS\\tFILTER';
        bcftools query -f '%QUAL\t%INFO/QD\t%INFO/SOR\t%INFO/FS\t%FILTER\n' ${params.time}.soft-filter.vcf.gz;
      } > filter_stats.tsv
  """
}


/*

Hard filter

This filtration is used to get 
a. biallelic SNPs
b. homozygous SNPs
*/

process HARD_FILTER {

  tag { "Perform hard filter" }

  label 'large'

  publishDir "${params.output}/variation/hard_filter", mode: 'copy'

  input:
      tuple path(vcf), path(vcf_index), path(contigs)

  output:
      tuple path("${params.time}.hard-filter.vcf.gz"), path("${params.time}.hard-filter.vcf.gz.tbi"), path("${params.time}.hard-filter.vcf.gz.csi"), emit: hard_filter_vcf
      path("${params.time}.hard-filter.stats.txt"), emit: hard_vcf_stats

  """
      # Create hard filter function
      function perform_hard_filter {

        bcftools view -m2 -M2 --trim-alt-alleles -O u --regions \${1} ${vcf} |\\
        bcftools filter -O u --set-GTs . --exclude 'FORMAT/FT ~ "DP_min_depth"' |\\
        bcftools filter -O u --exclude 'FILTER != "PASS"' |\\
        bcftools view -O v --min-af 0.000001 --max-af 0.999999 |\\
        bcftools view -O z --trim-alt-alleles > \${1}.vcf.gz        
      }

      export -f perform_hard_filter

      parallel --verbose perform_hard_filter {} < ${contigs}

      awk '{print \$0 ".vcf.gz"}' ${contigs} > contig_set.tsv
      bcftools concat -O z --file-list contig_set.tsv > ${params.time}.hard-filter.vcf.gz

      # Indexing vcf file
      bcftools index ${params.time}.hard-filter.vcf.gz
      bcftools index --tbi ${params.time}.hard-filter.vcf.gz
      bcftools stats -s- --verbose ${params.time}.hard-filter.vcf.gz > ${params.time}.hard-filter.stats.txt
  """
}







