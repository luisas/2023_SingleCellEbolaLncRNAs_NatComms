// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"

params.output_dir_name = "01_scRNA-Seq_inVivo_rhemac10"
params.primer = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
params.strandness = "FR"

params.output_dir_preliminary = "${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/"
params.output_dir = "${params.dirData}/02_scRNA-Seq_PBMCs/${params.output_dir_name}"


Channel.fromPath("${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")
       .into{ GtfChannel;GtfChannel2; GtfChannel3;GtfChannel4; GtfChannel5;   }


Channel.fromPath("${params.output_dir}/02_star/*/*/*/*/*/*merged.bam")
              .ifEmpty('bam files directory is empty')
              .map { tuple(it.baseName.split('_')[0],
                           it.baseName.split('_')[2],
                           it.baseName.split('_')[1],
                           it.baseName.split('_')[3].split("\\.")[0],
                           "std",
                           it.baseName.split('_')[0]+ "_" + it.baseName.split('_')[1] + "_" + it.baseName.split('_')[2] + "_" + it.baseName.split('_')[3],
                           it ) }.into{ bams_1; bams_2; bams_3;  }

bams_2.subscribe{ println "$it" }

params.scripts="${baseDir}/scripts/"

cell_selection_script = Channel.fromPath("${params.scripts}/cellselection.R").collect()
process TagReadWithGeneFunction{

  cpus 12

  storeDir "${params.output_dir}/03_DropSeqPreProcessing/$animal_id/$hpi/$exp/$replicate/$preprocessing"

  input:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file(bam) from bams_1
  file(gtf) from GtfChannel5.collect()

  output:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file("${complete_id}_gene_exon_tagged.bam") into bams_tagged


  script:
  """
  TagReadWithGeneFunction I=${bam} \
                              O=${complete_id}_gene_exon_tagged.bam \
                              ANNOTATIONS_FILE=${gtf} \
                              USE_STRAND_INFO=true
  """

}

process DetectBeadSubstitutionError{

  cpus 12

  storeDir "${params.output_dir}/03_DropSeqPreProcessing/$animal_id/$hpi/$exp/$replicate/$preprocessing"

  input:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file(bam) from bams_tagged

  output:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file("${complete_id}_clean.bam") into bams_tagged_clean
  file "*" into stats_beadsubstitution

  script:
  """
  DetectBeadSynthesisErrors \
      I=${bam} \
      O=${complete_id}_clean.bam \
      REPORT=${complete_id}.indel_report.txt \
      OUTPUT_STATS=${complete_id}.synthesis_stats.txt \
      SUMMARY=${complete_id}.synthesis_stats.summary.txt \
      PRIMER_SEQUENCE=${params.primer}
  """
}

process BamTagHistogram{

  cpus 12

  storeDir "${params.output_dir}/04_CellSelection/$animal_id/$hpi/$exp/$replicate/$preprocessing"

  input:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file(bam) from bams_tagged_clean

  output:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file(bam), file("${complete_id}_out_cell_readcounts.txt.gz") into BamTagHistogramChannel


  script:
  """
  BamTagHistogram I=${bam} O=${complete_id}_out_cell_readcounts.txt.gz TAG=XC
  """
}


process CellSelection{

  cpus 4

  storeDir "${params.output_dir}/03_DropSeqPreProcessing/$animal_id/$hpi/$exp/$replicate/$preprocessing"
  input:
  file cell_selection_script
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file(bam), file(cell_readcounts) from BamTagHistogramChannel

  output:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file(bam), file("${complete_id}_barcodes.txt") into CellSelectionChannel

  script:
  """
  Rscript ${cell_selection_script} ${cell_readcounts} ${complete_id}_barcodes.txt
  """
}

CellSelectionChannel.into{CellSelectionChannel1; CellSelectionChannel2}

// OUTPUT_READS_INSTEAD=Boolean  Output number of reads instead of number of unique molecular barcodes.  Default value:
// false. This option can be set to 'null' to clear the default value. Possible values:{true, false}
process DigitalExpressionMatrix_UMI{

  cpus 24

  storeDir "${params.output_dir}/04_DigitalExpressionMatrix/$animal_id/$hpi/$exp/$replicate/$preprocessing"

  input:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file(bam), file(barcodes) from CellSelectionChannel1

  output:
  file "*" into DigitalExpressionUMIChannel

  script:
  """
  java -jar /Drop-seq_tools-2.3.0/jar/dropseq.jar DigitalExpression \
          I=${bam}\
          O=${complete_id}.dge.txt.gz\
          SUMMARY=${complete_id}.summary.txt\
          OUTPUT_LONG_FORMAT=${complete_id}.dge.summary.txt\
          OUTPUT_READS_INSTEAD=false\
          CELL_BC_FILE=${barcodes}
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
