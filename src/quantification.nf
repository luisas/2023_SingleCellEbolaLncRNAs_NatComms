#!/usr/bin/env nextflow

//params.output_dir = "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all"
params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_Ebola-RNASeq_all"




//params.dataset= "02_RNA-Seq_BatchZyagen"
//params.htseqsense="yes"
params.dataset = "02_RNA-Seq_external"
params.htseqsense="reverse"


// Fastq1 have no lanes
files = Channel.fromPath("${params.output_dir}/${params.dataset}/03_hisat/*/*/*/*/*.UMI.f3.q60.bam")
              .ifEmpty('fastq files directory is empty')
              .map { tuple(it.baseName.split('_')[0],
                              it.baseName.split('_')[2],
                              it.baseName.split('_')[1],
                              it.baseName.split('_')[3].split("\\.")[0],
                              it.baseName.split('_')[0]+"_"+
                                              it.baseName.split('_')[1]+"_"+
                                              it.baseName.split('_')[2]+"_"+
                                              it.baseName.split('_')[3].split("\\.")[0],
                             it) }

files.into{ filtered_merged_bams_4; fastq_files_for_mapping; printing}

printing.subscribe{ println it }





//params.gtf_ref_merged =  "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all/03_novel_lncrnas/00_gtf_filtering_step_01/rheMac10_EBOV-Kikwit_and_both_novel.gtf"
params.gtf_ref_merged =  "/gpfs/projects/bsc83/Data/Ebola/01_Ebola-RNASeq_all/03_novel_lncrnas/02_final_catalogue/rheMac10_EBOV-Kikivit_and_novelcatalogue_with_names.gtf"

Channel.fromPath("${params.gtf_ref_merged}").set{ gtfANDnovel }

process getHTseqCountALL{
  cpus 8
  storeDir "${params.output_dir}/04_Quantification_correct_with_names/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"

  input:
  file gtf from gtfANDnovel.collect()
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bam) from filtered_merged_bams_4

  output:
  set file("${complete_id}_namesorted.bam"), file("${bam_prefix}.HTseq.gene_counts.tab") into htseqCountsALL_channel

  script:
  bam_prefix = get_file_name_no_extension(bam.name)
  """
  samtools sort ${bam} -@ ${task.cpus} -O bam -o ${complete_id}_namesorted.bam
  # Calc the counts for the umi_dedup
  htseq-count -f bam -r name -s ${params.htseqsense} -t exon -i gene_id ${complete_id}_namesorted.bam ${gtf} > ${bam_prefix}.HTseq.gene_counts.tab
  """
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
