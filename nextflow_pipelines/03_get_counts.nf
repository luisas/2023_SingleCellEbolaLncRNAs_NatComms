#!/usr/bin/env nextflow

/*
*
* Nextflow pipeline for filtering and quantifying the novel lncRNAs.
*
*/
log.info "=============================================="
log.info " Quantification and Filtering  "
log.info "=============================================="

// ------------------------------------------------------------
// ------------ INPUT PARAMETERS ------------------------------
// ------------------------------------------------------------

params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/"
params.dataset= "02_RNA-Seq_ribodepl"
params.data_dir = "${params.output_dir}/${params.dataset}"
params.study="*"
params.strandness = "FR"
if( "${params.strandness}" == "FR" ){
  params.htseqsense="yes"
  params.stringtiestrandness="fr"
}
if( "${params.strandness}" == "RF" ){
  params.htseqsense="reverse"
  params.stringtiestrandness="rf"
}

reference = Channel.fromPath("${params.output_dir}/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")

bams = Channel.fromFilePairs("${params.output_dir}/${params.dataset}/03_hisat/*${params.study}/*/*/*/*.UMI.f3.q60.umi_dedup.{bam,bam.bai}")
                .ifEmpty("No bams found")
                .map { tuple(it[0].split('_')[0],
                             it[0].split('_')[2],
                             it[0].split('_')[1],
                             it[0].split('_')[3].split("\\.")[0],
                             it[0].split('_')[0]+"_"+
                                             it[0].split('_')[1]+"_"+
                                             it[0].split('_')[2]+"_"+
                                             it[0].split('_')[3].split("\\.")[0],
                             it[1] ) }


process htseqcounts{

  storeDir "${params.output_dir}/04b_quantification_htseq/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  file gtf from reference.collect()
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bampair) from bams

  output:
  file "${complete_id}.HTseq.gene_counts.tab" into htseqCountsMerged_channel

  script:
  """
  # Calc the counts for the umi_dedup
  htseq-count -f bam -r pos -s ${params.htseqsense} -t exon -i gene_id ${bampair[0]} ${gtf} > ${complete_id}.HTseq.gene_counts.tab
  """
}



 workflow.onComplete {
 	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
 }

 /*   -------------------------------
 *           Groovy Functions
 *    -------------------------------
 */

 def get_file_name_no_extension(String filename){
   return filename.split("\\.").init().join('.')
 }
