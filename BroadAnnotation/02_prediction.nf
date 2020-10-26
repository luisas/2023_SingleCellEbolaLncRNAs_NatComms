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

params.prefix_data = "/gpfs/projects/bsc83/Data"
params.output_dir = "${params.prefix_data}/Ebola/99_BroadAnnotation/"
params.output_dir_preliminary = "${params.prefix_data}/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/"
// Create 2 different batch files
//params.htseqsense="yes"
//params.umis = "true"
//params.strandrule="1++,1--,2+-,2-+"
params.prefix_label = "ribodepleted"

files = Channel.fromPath("/gpfs/projects/bsc83/Data/Ebola/00_RawData/BroadTranscriptomesComplete/bams/*.bam")
                .ifEmpty("No bams found")
                .map { tuple(it.baseName, it) }


files.into{ filtered_merged_bams_4; fastq_files_for_mapping; printing; filtered_merged_bams_5}
//printing.subscribe{ println it }


params.data_dir = "${params.output_dir}/"
params.lnc_novel_compared_to_ref= "${params.data_dir}/01_stringtie_assembly_merged/01_gffCompare/merged.annotated.gtf"
params.assembly ="${params.data_dir}/01_stringtie_assembly_merged/stringtie_merged_reference_guided.gtf"
params.reference_annotated = "${params.output_dir_preliminary}/gene_annotations/UCSC/rheMac10_EBOV-Kikwit_UCSC.gtf"
reference_assembly = Channel.fromPath("${params.output_dir_preliminary}/reference_assembly/rheMac10_EBOV-Kikwit_UCSC.fa")

// // SCRIPTS
extractNovel_script = Channel.fromPath("${baseDir}/scripts/03_concordant_and_merge.R").collect()
prefilter_script = Channel.fromPath("${baseDir}/scripts/01_prefilter_script.R").collect()
expression_script = Channel.fromPath("${baseDir}/scripts/02_filterExpressed.R").collect()
lnc_novel_compared_to_ref = Channel.fromPath("${params.lnc_novel_compared_to_ref}").collect()
mrnas_predicted = Channel.fromPath("${params.mrnas_predicted}").collect()
assembly = Channel.fromPath("${params.assembly}").collect()
gtf2bed = Channel.fromPath("${baseDir}/scripts/gtf2bed.R").collect()
gtf_ref_1 = Channel.fromPath("${params.reference_annotated}").collect()
ref = Channel.fromPath("${params.reference_annotated}").collect()

cpatmodel= Channel.fromPath("${params.output_dir_preliminary}/Human_logitModel.RData").collect()
cpathexamer= Channel.fromPath("${params.output_dir_preliminary}/Human_Hexamer.tsv").collect()
cpat_files = Channel.fromPath("${params.output_dir_preliminary}/cpat*").collect()

gtf_ref_1.into{gtf_ref; gtf_ref_2; }
// /*
// * Extract only the novel (unkown and antisense) lncRNAs compared to the reference.
// * Also remove non concordant predictions.
// */
output_dir_sub="03_novel_lncRNAs_list"

process prefilter{

   storeDir "${params.output_dir}/${output_dir_sub}/00_prefilter_candidates"

   input:
   file prefilter_script
   file assembly
   file lnc_novel_compared_to_ref
   file ref

   output:
   file("prefilter_candidates.gtf") into (prefiter_candidates_channel, prefiter_candidates_channel2)

   script:
   """
   Rscript ${prefilter_script} ${lnc_novel_compared_to_ref} ${assembly} ${ref} prefilter_candidates.gtf
   """
}
//
// /*
// * Append novel concordant to the reference before quantification
// */
process mergeWithAnnotation {
     storeDir "${params.output_dir}/${output_dir_sub}/00_prefilter_candidates/"

     input:
     file gtf_ref
     file novel from prefiter_candidates_channel

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

   storeDir "${params.output_dir}/${output_dir_sub}/01_quantification_for_filtering/stringtie"

   input:
   file gtf from novel_and_reference.collect()
   set complete_id,
       file(bampair) from filtered_merged_bams_4

   output:
   file("${complete_id}.gene_abundances.tsv") into fpkm_channel

   script:
   """
   # Calc the counts for the umi_dedup
   stringtie -eB -G ${gtf} ${bampair[0]} -A ${complete_id}.gene_abundances.tsv
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
   storeDir "${params.output_dir}/${output_dir_sub}/02_novel_expressed/"

   input:
   file expression_script
   file prefilter_candidates from prefiter_candidates_channel2
   file expression

   output:
   file("candidates.gtf") into novel_expressed

   script:
   """
   mkdir expression_dir
   mv ${expression} expression_dir
   Rscript ${expression_script} expression_dir ${prefilter_candidates} candidates.gtf
   """
 }


novel_expressed.into{novel_expressed1; novel_expressed2; novel_expressed3; }


process gtf2bed{
   storeDir "${params.output_dir}/${output_dir_sub}/"

   input:
   file gtf2bed
   file candidates from novel_expressed1

   output:
   file("${candidates.baseName}.bed12") into candidatesbed

   script:
   """
   Rscript ${gtf2bed} ${candidates} ${candidates.baseName}.bed12
   """
}

process getFasta{
   storeDir "${params.output_dir}/${output_dir_sub}/"

   input:
   file gtf2bed
   file reference_assembly
   file candidates from candidatesbed

   output:
   file("${candidates.baseName}.fa") into candidatesfa

   script:
   """
   bedtools getfasta -fi ${reference_assembly} -bed ${candidates} -name -fo ${candidates.baseName}.fa -s
   """
}


candidatesfa.into{candidatesfa1; candidatesfa2; candidatesfa3; candidatesfa4; }

process CPC2{
  storeDir "${params.output_dir}/${output_dir_sub}/03_predictions/CPC2"
   input:
  file candidates from candidatesfa1
   output:
   file("*") into cpc2_all
   file("cpc2_pred.txt") into cpc2

   script:
   """
   python2 \${CPC_HOME}/CPC2.py -i ${candidates} -o cpc2_pred
   """
}


process CPAT{
   storeDir "${params.output_dir}/${output_dir_sub}/03_predictions/CPAT"

   input:
   file candidates from candidatesfa2
   file cpatmodel
   file cpathexamer
   file cpat_files

   output:
   file("*") into cpat_all
   file("cpat_pred.ORF_prob.tsv") into cpat

   script:
   """
   cpat.py -g ${candidates} -o cpat_pred  -d ${cpatmodel} -x ${cpathexamer}
   """
}



process CNIT{
    storeDir "${params.output_dir}/${output_dir_sub}/03_predictions/CNIT"

    input:
    file candidates from candidatesfa3

    output:
    file("*") into cnit_all
    file("cnit_pred/CNCI2.index") into cnit
    script:
    """
    python \${CNIT_HOME}/CNCI2.py -f ${candidates} -o cnit_pred -m "ve"
    """
}


process getFinalAnnotation {
     storeDir "${params.output_dir}/03_novel_lncRNAs_list/"

     input:
     file cpatpred from cpat
     file cnitpred from cnit
     file cpc2pred from cpc2
     file gtf_ref_2
     file extractNovel_script
     file candidates from novel_expressed
     output:
     set file("rheMac10_EBOV_and_novel.gtf"), file("rheMac10_EBOV_and_novel_genenames.gtf") into (novel_final)

     script:
     """
     Rscript ${extractNovel_script} ${cpc2pred} ${cpatpred} ${cnitpred} ${candidates} ${gtf_ref_2} rheMac10_EBOV_and_novel.gtf rheMac10_EBOV_and_novel_genenames.gtf
     """
 }



process stringtie{
   cpus 1
   storeDir "${params.output_dir}/04_quantification/"

   input:
   set file(gtf), file(gtf_names) from novel_final.collect()
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
