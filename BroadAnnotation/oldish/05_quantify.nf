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

params.prefix_data = "/gpfs/projects/bsc83/Data/Ebola"
params.output_dir_name = "99_BroadAnnotation_rerunningassemly_old"
params.output_dir = "${params.prefix_data}/${params.output_dir_name}/"
// Create 2 different batch files
//params.htseqsense="yes"
//params.umis = "true"
//params.strandrule="1++,1--,2+-,2-+"


files = Channel.fromPath("${params.prefix_data}/00_RawData/BroadTranscriptomesComplete/bams/*.bam")
                .ifEmpty("No bams found")
                .map { tuple(it.baseName, it) }


files.into{ bams; filtered_merged_bams_4; fastq_files_for_mapping; printing; filtered_merged_bams_5}
printing.subscribe{ println it }


// /*
// * Extract only the novel (unkown and antisense) lncRNAs compared to the reference.
// * Also remove non concordant predictions.
// */

params.output_dir_sub = "03_novel_lncRNAs_list"
output_dir_sub= "${params.output_dir_sub}"


novel_final = Channel.fromPath("${params.output_dir}/${params.output_dir_sub}/rheMac10_novel.gtf")

process stringtie{
   cpus 1
   storeDir "${params.output_dir}/${output_dir_sub}/04_quantification_new_new/$complete_id"

   input:
   set file(gtf_names) from novel_final.collect()
   set complete_id,
       file(bampair) from filtered_merged_bams_5

   output:
   file("*") into expression_channel

   script:
   """
   # Calc the counts for the umi_dedup
   stringtie -eB -G ${gtf_names} ${bampair[0]} -A ${complete_id}.gene_abundances.tsv
   mv t_data.ctab ${complete_id}_t_data.ctab
   """
}



 workflow.onComplete {
 	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
 }

 workflow.onComplete {
 	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
 }

 /*   -------------------------------
 *           Groovy Functions
 *    -------------------------------
 */

 def remove_lane_from_id(String id){
   return id.split("_").init().join("_")
 }

 def get_file_name_no_extension(String filename){
   return filename.split("\\.").init().join('.')
 }
