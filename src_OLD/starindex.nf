// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"


//params.dirData = "/home/luisas/Desktop/cluster/data/"
//params.dirProj = "/home/luisas/Desktop/cluster/proj/Ebola"
params.output_dir_preliminary = "${params.dirData}/01_Ebola-RNASeq_all/01_PreliminaryFiles_rheMac8/"


Channel.fromPath("${params.output_dir_preliminary}/gene_annotations/rheMac8_EBOV-Kikiwit.gtf")
       .into{ GtfChannel;   }
Channel.fromPath("${params.output_dir_preliminary}/reference_assembly/rheMac8_EBOV-Kikiwit.fa")
       .into{ ReferenceChannel;}






log.info "=============================================="
log.info " Star Index generatiomn   "
log.info "=============================================="



process create_star_indexes{

  storeDir "${params.output_dir_preliminary}/indexes/star"
  cpus 48
  label "big_mem"

  input:
  file assembly from ReferenceChannel
  file annotation from GtfChannel

  output:
  file "star" into star_indexes

  script:
  """
  mkdir star
  STAR --runMode genomeGenerate \\
    --sjdbGTFfile $annotation \\
    --genomeDir star/ \\
    --genomeFastaFiles $assembly \\
    --runThreadN ${task.cpus}
  """
}


workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}
