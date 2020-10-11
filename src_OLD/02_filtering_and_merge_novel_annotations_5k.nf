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

// Create 2 different batch files
params.dataset= "02_RNA-Seq_ribodepl"
//params.htseqsense="yes"
//params.umis = "true"
//params.strandrule="1++,1--,2+-,2-+"
params.prefix_label = "ribodepleted"

//params.dataset = "02_RNA-Seq_external"
//params.htseqsense="reverse"
//params.umis = "false"
//params.strandrule="1+-,1-+,2++,2--"
//params.prefix_label = "polya"



files = Channel.fromFilePairs("${params.output_dir}/${params.dataset}/03_hisat/*/*/*/*/*.UMI.f3.q60.umi_dedup.{bam,bam.bai}")
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


files.into{ filtered_merged_bams_4; fastq_files_for_mapping; printing; filtered_merged_bams_5}
printing.subscribe{ println it }

params.data_dir = "${params.output_dir}/${params.dataset}"
params.lnc_novel_compared_to_ref= "${params.data_dir}/05_feelNC_5k/01_gffcompare/merged.annotated.gtf"
params.mrnas_predicted ="${params.data_dir}/05_feelNC_5k/feelnc_codpot_out/candidate_lncRNA.gtf.mRNA.gtf"
params.lnc_predicted ="${params.data_dir}/05_feelNC_5k/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf"
params.reference_annotated = "${params.output_dir}/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit.gtf"


// // SCRIPTS
extractNovel_script = Channel.fromPath("${baseDir}/scripts/03_addNovel.R").collect()
novel_concordant_script = Channel.fromPath("${baseDir}/scripts/01_filter_novel_concordant.R").collect()
expression_script = Channel.fromPath("${baseDir}/scripts/02_filterExpressed.R").collect()
lnc_novel_compared_to_ref = Channel.fromPath("${params.lnc_novel_compared_to_ref}").collect()
mrnas_predicted = Channel.fromPath("${params.mrnas_predicted}").collect()
lnc_predicted = Channel.fromPath("${params.lnc_predicted}").collect()
gtf_ref_1 = Channel.fromPath("${params.reference_annotated}").collect()
ref = Channel.fromPath("${params.reference_annotated}").collect()
gtf_ref_1.into{gtf_ref; gtf_ref_2; }
// /*
// * Extract only the novel (unkown and antisense) lncRNAs compared to the reference.
// * Also remove non concordant predictions.
// */


process novelANDconcordant{

   storeDir "${params.output_dir}/03_novel_lncRNAs_list_5k/00_all_novels/"

   input:
   file novel_concordant_script
   file lnc_predicted
   file lnc_novel_compared_to_ref
   file mrnas_predicted
   file ref

   output:
   file("novel_rhemac10_concordant_${params.prefix_label}.gtf") into (novel_concordant_channel, novel_concordant_channel2)

   script:
   """
   Rscript ${novel_concordant_script} ${lnc_novel_compared_to_ref} ${lnc_predicted} ${mrnas_predicted} ${ref} novel_rhemac10_concordant_${params.prefix_label}.gtf
   """
}
//
// /*
// * Append novel concordant to the reference before quantification
// */
process mergeWithAnnotation {
     storeDir "${params.output_dir}/03_novel_lncRNAs_list_5k/00_all_novels/"

     input:
     file gtf_ref
     file novel from novel_concordant_channel

     output:
     file("ref_novel_${params.prefix_label}.gtf") into (novel_and_reference, novel_and_reference2)

     script:
     """
     cat ${gtf_ref} > ref_novel_${params.prefix_label}.gtf
     cat ${novel} >> ref_novel_${params.prefix_label}.gtf
     """
 }
//
// /*
// * Quantify Expression with stringtie.
// */
 process expression_stringtie{

   storeDir "${params.output_dir}/03_novel_lncRNAs_list_5k/01_quantification_for_filtering/$dataset_name/$tissue/$dayPostInfection/$sample/stringtie"

   input:
   file gtf from novel_and_reference.collect()
   set dataset_name, dayPostInfection, tissue, sample, complete_id,
       file(bampair) from filtered_merged_bams_4

   output:
   file("${complete_id}.gene_abundances.tsv") into fpkm_channel

   script:
   """
   # Calc the counts for the umi_dedup
   stringtie -eB  --fr -G ${gtf} ${bampair[0]} -A ${complete_id}.gene_abundances.tsv
   mv t_data.ctab ${complete_id}_t_data.ctab
   """
 }
//
//
 expression = fpkm_channel.collect()
 expression_script = Channel.fromPath("${baseDir}/scripts/02_filterExpressed.R").collect()



 /*
 * Quantify Expression with stringtie.
 */
 process filterExpression{
   storeDir "${params.output_dir}/03_novel_lncRNAs_list_5k/02_novel_expressed/"

   input:
   file expression_script
   file novelconcordant from novel_concordant_channel2
   file expression

   output:
   file("novel_expressed_${params.prefix_label}.gtf") into novel_expressed

   script:
   """
   mkdir expression_dir
   mv ${expression} expression_dir
   Rscript ${expression_script} expression_dir ${novelconcordant} novel_expressed_${params.prefix_label}.gtf
   """
 }





process getNovel {
    storeDir "${params.output_dir}/03_novel_lncRNAs_list_5k/"

    input:
    file gtf_ref_2
    file extractNovel_script
    file novel_ribodepl from novel_expressed

    output:
    set file("rheMac10_EBOV_and_novel.gtf"), file("rheMac10_EBOV_and_novel_genenames.gtf") into (novel_final)

    script:
    """
    Rscript ${extractNovel_script} ${novel_ribodepl} ${gtf_ref_2} rheMac10_EBOV_and_novel.gtf rheMac10_EBOV_and_novel_genenames.gtf
    """
}



process stringtie{
  cpus 1
  storeDir "${params.output_dir}/04_quantification/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set file(gtf), file(gtf_names) from novel_final.collect()
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bampair) from filtered_merged_bams_5

  output:
  file("*") into expression_channel

  script:
  """
  # Calc the counts for the umi_dedup
  stringtie -eB  --fr -G ${gtf_names} ${bampair[0]} -A ${complete_id}.gene_abundances.tsv
  mv t_data.ctab ${complete_id}_t_data.ctab
  """
}


//
//
 workflow.onComplete {
 	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
 }
//
// //  #FPKM_count.py -i ${bampair[0]} -r ${bed12} -o ${complete_id} --strand ${params.strandrule} --only-exonic
// /*   -------------------------------
// *           Groovy Functions
// *    -------------------------------
// */
//
 def remove_lane_from_id(String id){
   return id.split("_").init().join("_")
 }
//
 def get_file_name_no_extension(String filename){
   return filename.split("\\.").init().join('.')
 }
