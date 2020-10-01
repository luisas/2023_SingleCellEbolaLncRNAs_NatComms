



expression_script = Channel.fromPath("${baseDir}/scripts/02_filterExpressed.R").collect()

expression = Channel.fromPath("${params.output_dir}/${params.dataset}/03_hisat/*/*/*/*/*.UMI.f3.q60.{bam,bam.bai}")
                .ifEmpty("No bams found")


params.prefix_label = "ribodepleted"
novel_concordant_channel2 = Channel.fromPath("${params.output_dir}/03_novel_lncRNAs_list/00_all_novels/novel_rhemac10_concordant_${params.prefix_label}.gtf")


expression.into{ expression1; expression2; }

expression1.subscribe{ println "$it"}


/*
* Quantify Expression with stringtie.
*/
process filterExpression{
  storeDir "${params.output_dir}/03_novel_lncRNAs_list/02_novel_expressed/"

  input:
  file expression_script
  file novelconcordant from novel_concordant_channel2
  file expression2

  output:
  file("novel_expressed_${params.prefix_label}.gtf") into novel_expressed

  script:
  """
  mkdir expression_dir
  mv ${expression} expression_dir
  Rscript ${expression_script} expression_dir ${novelconcordant} novel_expressed_${params.prefix_label}.gtf
  """
}
