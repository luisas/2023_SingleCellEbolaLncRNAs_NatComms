// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"

//params.dirData = "/home/luisas/Desktop/cluster/data/"
//params.dirProj = "/home/luisas/Desktop/cluster/proj/Ebola"

params.dataset_bam_dir = "${params.dirData}/00_RawData/scRNAseq_exvivo_alin_bams"
params.output_dir = "${params.dirData}/01_Ebola-RNASeq/03_scRNA-Seq"


Channel.fromPath("${params.dataset_bam_dir}/*.bam")
              .ifEmpty('bam files directory is empty')
              .map { tuple(it.baseName.split('\\.')[0],
                           it.baseName.split('\\.')[1],
                           it.baseName.split('\\.')[3],
                           it.baseName.split('\\.')[6],
                           it.baseName.split('\\.')[0]+ "_" + it.baseName.split('\\.')[6] + "_" + it.baseName.split('\\.')[1] + "_" + it.baseName.split('\\.')[3],
                           it ) }.into{ bams_1; bams_2; bams_3;  }

Channel.fromPath("${params.dirData}/01_Ebola-RNASeq/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV-Kikwit_and_novel_with_names.gtf")
       .into{ GtfChannel; }

cell_selection_script = Channel.fromPath("${baseDir}/scripts/cellselection.R").collect()
// We need to estimate how manty cells we want to extract


process TagReadWithGeneFunction{

  storeDir "${params.output_dir}/00_Tag/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id, file(bam) from bams_1
  file(gtf) from GtfChannel.collect()

  output:
  set animal_id, hpi, exp, replicate, complete_id, file("${complete_id}_gene_exon_tagged.bam") into bams_tagged


  script:
  """
  TagReadWithGeneFunction I=${bam} \
                              O=${complete_id}_gene_exon_tagged.bam \
                              ANNOTATIONS_FILE=${gtf}
  """

}


process BamTagHistogram{

  storeDir "${params.output_dir}/00_BamTagHistogram/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id, file(bam) from bams_tagged

  output:
  set animal_id, hpi, exp, replicate, complete_id, file(bam), file("${complete_id}_out_cell_readcounts.txt.gz") into BamTagHistogramChannel


  script:
  """
  BamTagHistogram I=${bam} O=${complete_id}_out_cell_readcounts.txt.gz TAG=XC
  """
}

process CellSelection{

  storeDir "${params.output_dir}/00_BamTagHistogram/$animal_id/$hpi/$exp/$replicate"
  input:
  file cell_selection_script
  set animal_id, hpi, exp, replicate, complete_id, file(bam), file(cell_readcounts) from BamTagHistogramChannel

  output:
  set animal_id, hpi, exp, replicate, complete_id, file(bam), file("${complete_id}_barcodes.txt") into CellSelectionChannel

  script:
  """
  Rscript ${cell_selection_script} ${cell_readcounts} ${complete_id}_barcodes.txt
  """
}


CellSelectionChannel.into{CellSelectionChannel1; CellSelectionChannel2}


// // OUTPUT_READS_INSTEAD=Boolean  Output number of reads instead of number of unique molecular barcodes.  Default value:
// // false. This option can be set to 'null' to clear the default value. Possible values:{true, false}
process DigitalExpressionMatrix_UMI{

  cpus 24
  storeDir "${params.output_dir}/01_digital_expression/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id, file(bam), file(barcodes) from CellSelectionChannel1

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


//
// process DigitalExpressionMatrix_READS{
//
//   storeDir "${params.output_dir}/01_digital_expression/$animal_id/$hpi/$exp/$replicate"
//
//   input:
//   set animal_id, hpi, exp, replicate, complete_id, file(bam), file(barcodes) from CellSelectionChannel2
//
//   output:
//   file "*" into DigitalExpressionREADSChannel
//
//   script:
//   """
//   java -jar /Drop-seq_tools-2.3.0/jar/dropseq.jar DigitalExpression \
//           I=${bam}\
//           O=${complete_id}_reads_.dge.txt.gz\
//           SUMMARY=${complete_id}_reads_.summary.txt\
//           OUTPUT_LONG_FORMAT=${complete_id}_reads_.dge.summary.txt\
//           OUTPUT_READS_INSTEAD=true\
//           CELL_BC_FILE=${barcodes}
//   """
// }



// Convert matrix
