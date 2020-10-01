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

params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation"


params.data_dir = "${params.output_dir}/${params.dataset}"
params.lnc_novel_compared_to_ref= "${params.data_dir}/05_feelNC_prediction/feelnc_gencode_linc/01_gffcompare/merged.annotated.gtf"
params.mrnas_predicted ="${params.data_dir}/05_feelNC_prediction/feelnc_gencode_linc/feelnc_codpot_out/candidate_lncRNA.gtf.mRNA.gtf"
params.reference_annotated = "${params.output_dir}/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit.gtf"

extractNovelFromPolyA_script = Channel.fromPath("${baseDir}/scripts/03_getNovelBetweenPredictions.R").collect()
novel_ribodepl = Channel.fromPath("${params.output_dir}/03_novel_lncRNAs_list/02_novel_expressed/novel_expressed_ribodepleted.gtf").collect()
novel_polyA = Channel.fromPath("${params.output_dir}/03_novel_lncRNAs_list/02_novel_expressed/novel_expressed_polya.gtf").collect()
gtf_ref = Channel.fromPath("${params.output_dir}/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit.gtf").collect()

/*
* Extract only the novel (unkown and antisense) lncRNAs compared to the reference.
* Also remove non concordant predictions.
*/
process gffCompare{

  storeDir "${params.output_dir}/03_novel_lncRNAs_list/03_polyA_vs_ribodepl/"

  input:
  file novel_polyA
  file novel_ribodepl

  output:
  file("ribodeplVSpolyA.annotated.gtf") into gffcomparech

  script:
  """
  gffcompare -R -r ${novel_polyA} -o ribodeplVSpolyA ${novel_ribodepl}
  """
}

/*
* Append novel concordant to the reference before quantification
*/
process getNovel_fromPolyA {
    storeDir "${params.output_dir}/03_novel_lncRNAs_list/"

    input:
    file gtf_ref
    file ribodeplVSpolyA from gffcomparech
    file extractNovelFromPolyA_script
    file novel_polyA
    file novel_ribodepl

    output:
    set file("rheMac10_EBOV_and_novel.gtf"), file("rheMac10_EBOV_and_novel_genenames.gtf") into (novel_and_reference, novel_and_reference2)

    script:
    """
    Rscript ${extractNovelFromPolyA_script} ${ribodeplVSpolyA} ${novel_ribodepl} ${novel_polyA} ${gtf_ref} rheMac10_EBOV_and_novel.gtf rheMac10_EBOV_and_novel_genenames.gtf
    """
}


workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

//  #FPKM_count.py -i ${bampair[0]} -r ${bed12} -o ${complete_id} --strand ${params.strandrule} --only-exonic
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
