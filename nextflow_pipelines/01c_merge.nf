#!/usr/bin/env nextflow

//params.output_dir = "/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all"
params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_ribodepl"

// Reference GTF
params.reference_annotated = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit.gtf"
gtfChannel = Channel.fromPath("${params.reference_annotated}")
gtfChannel.into{gtfChannel1; gtfChannel2 }

// Stringtie Files
stringTie_channel = Channel.fromPath("${params.output_dir}/04_stringtie/*/*/*/*/*stringtie.gtf")


process StringTie_Merge_Reference_Guided{

  cpus 48
  storeDir "${params.output_dir}/04_stringtie/"

  input:
  file(stringtie_gtfs) from stringTie_channel.collect()
  file reference_gtf from gtfChannel1

  output:
  file "stringtie_merged_reference_guided.gtf" into (merged_denovo_assmebly, merged_de_novo_assembly_2)

  script:
  """
    stringtie --merge -p ${task.cpus} -o stringtie_merged_reference_guided.gtf -G ${reference_gtf} ${stringtie_gtfs}
  """
}


/*
* Compare stringtie output with reference GTF
*/
process gffCompare2{

  storeDir "${params.output_dir}/04_stringtie_gffcompare"

  input:
  file merged_gtf from merged_denovo_assmebly
  file reference_gtf from gtfChannel2

  output:
  file("merged*") into gff_compare_output_channel2

  script:
  """
  gffcompare -R -r ${reference_gtf} -o merged ${merged_gtf}
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
