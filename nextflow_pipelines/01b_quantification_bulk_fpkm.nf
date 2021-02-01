#!/usr/bin/env nextflow

//params.output_dir = "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all"
params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation"



params.dataset= "02_RNA-Seq_BatchZyagen"
params.htseqsense="yes"
params.umis = "true"
params.strandrule="1++,1--,2+-,2-+"
params.prefix_label = "ribodepleted"

//params.dataset = "02_RNA-Seq_external"
//params.htseqsense="reverse"
//params.umis = "false"
//params.strandrule="1+-,1-+,2++,2--"
//params.prefix_label = "polya"


if("${params.umis}" == "false"){
  files = Channel.fromFilePairs("${params.output_dir}/${params.dataset}/03_hisat/*/*/*/*/*.UMI.f3.q60.{bam,bam.bai}")
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
}else{
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
}

files.into{ filtered_merged_bams_4; fastq_files_for_mapping; printing}

printing.subscribe{ println it }


params.gtf= "/gpfs/projects/bsc83/Data/gene_annotation/ensembl_release100/rheMac10/Macaca_mulatta.Mmul_10.100.gtf"


process stringtie{
  cpus 1
  storeDir "${params.output_dir}/04_quantification/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  file gtf from gtf.collect()
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bampair) from filtered_merged_bams_4

  output:
  file("*") into fpkm_channel

  script:
  """
  # Calc the counts for the umi_dedup
  stringtie -eB  --fr -G ${gtf} ${bampair[0]} -A ${complete_id}.gene_abundances.tsv
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
