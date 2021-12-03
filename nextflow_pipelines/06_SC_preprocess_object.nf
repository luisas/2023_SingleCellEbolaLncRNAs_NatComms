// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"

params.output_dir_name = "01_scRNA-Seq_inVivo_rhemac10"
params.input_dir_name = "01_scRNA-Seq_inVivo_rhemac10"
params.output_dir = "${params.dirData}/02_scRNA-Seq_PBMCs/${params.output_dir_name}"
params.input_dir = "${params.dirData}/02_scRNA-Seq_PBMCs/${params.input_dir_name}"
params.output_prefix = "immune.combined"
params.scripts="${baseDir}/scripts/"
params.qc = "01_QC.R"
qc_script = Channel.fromPath("${params.scripts}/SC/${params.qc}").collect()
scrublet_script = Channel.fromPath("${params.scripts}/SC/02_scrublet_script.py").collect()

Channel.fromPath("${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")
       .into{ GtfChannel;GtfChannel2; GtfChannel3;GtfChannel4; GtfChannel5;   }

quantification_channel_1 = Channel.fromPath("${params.input_dir}/04_DigitalExpressionMatrix_OK/*/*/*/*/*/*.dge.txt.gz").collect()
quantification_channel_1.into{quantification_channel; quantification_channel_2}

quantification_channel_2.subscribe{ println "$it" }
println("${params.output_dir}")

process QC{

  cpus 1
  storeDir "${params.output_dir}/05_RObjects/01_QC"

  input:
  file quantification_channel
  file qc_script
  file ref from GtfChannel.collect()

  output:
  set file("${params.output_prefix}_qc.rds"), file("${params.output_prefix}_qc.mtx") into qc_filtered


  script:
  """
  mkdir data_dir
  mv ${quantification_channel} data_dir
  Rscript ${qc_script} data_dir ${ref} ${params.output_prefix}_qc.rds ${params.output_prefix}_qc.mtx
  """

}


process DoubletDetection{
  cpus 1
  storeDir "${params.output_dir}/05_RObjects/02_DoubletDetection"

  input:
  set file(rds), file(mtx) from qc_filtered
  file scrublet_script

  output:
  set file(rds), file(mtx) into qc_filtered_1
  set file("${params.output_prefix}_scrublet_mask.txt"), file("${params.output_prefix}_scrublet_score.txt") into doublet_removed

  script:
  """
  python3 ${scrublet_script} ${mtx} ${params.output_prefix}_scrublet_mask.txt ${params.output_prefix}_scrublet_score.txt
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
