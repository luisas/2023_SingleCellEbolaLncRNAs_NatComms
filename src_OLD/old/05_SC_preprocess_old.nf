// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"



params.output_dir_preliminary = "${params.dirData}/01_Ebola-RNASeq_all/01_PreliminaryFiles_rheMac10/"

params.dataset_bam_dir = "${params.dirData}/00_RawData/scRNAseq_exvivo_alin_bams"
params.output_dir = "${params.dirData}/01_Ebola-RNASeq_all/03_scRNA-Seq_rhemac8_mt"

//Channel.fromPath("${params.dataset_bam_dir}/EV0003.H024.fresh.a2.std.EX2.live-moi1e-1.bam")
Channel.fromPath("${params.dataset_bam_dir}/*.bam")
              .ifEmpty('bam files directory is empty')
              .map { tuple(it.baseName.split('\\.')[0],
                           it.baseName.split('\\.')[1],
                           it.baseName.split('\\.')[3],
                           it.baseName.split('\\.')[6],
                           it.baseName.split('\\.')[0]+ "_" + it.baseName.split('\\.')[6] + "_" + it.baseName.split('\\.')[1] + "_" + it.baseName.split('\\.')[3],
                           it ) }.into{ bams_1; bams_2; bams_3;  }

params.scripts="${baseDir}/scripts/"
gtfToGenePred_script_ch = Channel
                      .fromPath("${params.scripts}/gtfToGenePred")
genePredToBed_script_ch = Channel
                      .fromPath("${params.scripts}/genePredToBed")

//Channel.fromPath("${params.dirData}/01_Ebola-RNASeq_all/03_novel_lncrnas/02_final_catalogue/rheMac10_EBOV-Kikivit_and_novelcatalogue_with_names.gtf")
//       .into{ GtfChannel;GtfChannel2; GtfChannel3;GtfChannel4; GtfChannel5;   }
Channel.fromPath("${params.dirData}/01_Ebola-RNASeq_all/03_novel_lncrnas/02_final_catalogue/03_liftover/rheMac8_EBOV-Kikivit_and_novelcatalogue_with_names_ensembl.gtf")
       .into{ GtfChannel;GtfChannel2; GtfChannel3;GtfChannel4; GtfChannel5;   }

Channel.fromPath("${params.output_dir_preliminary}/gene_annotations/rheMac10_EBOV-Kikwit.gtf")
        .into{ GtfNotNovel;   }

Channel.fromPath("${params.output_dir_preliminary}/reference_assembly/rheMac10_EBOV-Kikwit.fa")
       .into{ ReferenceChannel;ReferenceChannel2;   }

Channel.fromPath("${params.dirData}/01_Ebola-RNASeq_all/01_PreliminaryFiles_rheMac10/reference_assembly/rheMac10_EBOV-Kikwit.dict")
       .into{ DictChannel; DictChannel2;DictChannel3; DictChannel4;  }

params.strandness = "FR"



cell_selection_script = Channel.fromPath("${baseDir}/scripts/cellselection.R").collect()
// We need to estimate how manty cells we want to extract

bams_2.subscribe{ println "$it" }

log.info "=============================================="
log.info " Metadata Generation  "
log.info "=============================================="

process ConvertToRefFlat{
  storeDir "${params.output_dir_preliminary}/gene_annotations"

  input:
  file(dict) from DictChannel
  file(gtf) from GtfChannel

  output:
  file("${prefix}.refFlat") into RefFlatChannel

  script:
  prefix = get_file_name_no_extension(gtf.name)
  """
  ConvertToRefFlat ANNOTATIONS_FILE=${gtf} \
                   SEQUENCE_DICTIONARY=${dict} \
                   OUTPUT=${prefix}.refFlat
  """
}

process ReduceGtf{
  storeDir "${params.output_dir_preliminary}/gene_annotations"

  input:
  file(dict) from DictChannel2
  file(gtf) from GtfChannel2

  output:
  file("${prefix}.reduced.gtf") into ReduceGTFChannel
  script:
  prefix = get_file_name_no_extension(gtf.name)
  """
  ReduceGtf SEQUENCE_DICTIONARY=${dict} \
            GTF=${gtf} \
            OUTPUT=${prefix}.reduced.gtf
  """
}

//need to add the other part
process convert_gtf_to_bed12{
  storeDir "${params.output_dir_preliminary}/gene_annotations"

  input:
  file annotation from GtfChannel4
  file gtfToGenePred_script from gtfToGenePred_script_ch
  file genePredToBed_script from genePredToBed_script_ch

  output:
  set file("${prefix}.bed"), file("${prefix}.bed12") into bed_channel

  script:
  prefix = get_file_name_no_extension(annotation.name)
  """
  ./${gtfToGenePred_script} ${annotation} ${prefix}.bed
  ./${genePredToBed_script} ${prefix}.bed ${prefix}.bed12
  """
}


process CreateIntervalFiles{

  storeDir "${params.output_dir_preliminary}/intervals"

  input:
  file(dict) from DictChannel3
  file(reducedgtf) from ReduceGTFChannel

  output:
  file("*") into IntervalsChannel

  script:
  prefix = get_file_name_no_extension(reducedgtf.name)
  """
  CreateIntervalsFiles SEQUENCE_DICTIONARY=${dict} \
                       REDUCED_GTF=${reducedgtf} \
                       PREFIX=${prefix} \
                       OUTPUT="."
  """
}


process create_star_indexes{

  storeDir "${params.output_dir_preliminary}/indexes/"
  cpus 16
  label "big_mem"

  input:
  file assembly from ReferenceChannel
  file annotation from GtfNotNovel

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

// Not doing preprocessing of bams cause they have already been preprocessed by Dylan
// preprocess is before mapping so independent from Annotation

log.info "=============================================="
log.info " DropSeq Analysis Pipeline  "
log.info "=============================================="

process BamToFastq{

  cpus 8
  tag "${complete_id}"
  storeDir "${params.output_dir}/00_fastq/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id, file(bam) from bams_1

  output:
  set animal_id, hpi, exp, replicate, complete_id, \
      file("${complete_id}.fastq.gz") into (fastqs_qc, fastqs_mapping)

  script:
  """
  picard-tools SamToFastq I=${bam} \
                          FASTQ=${complete_id}.fastq.gz

  """
}



/*
* Generate quality assessment with fastqc.
*/
process generate_fastqc{

  cpus 8
  label "rnaseq"
  tag "${complete_id}"
  storeDir "${params.output_dir}/01_fastqc/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id,
      file(fastq) from fastqs_qc

  output:
  file "*" into fastqcs

  script:
  """
  fastqc ${fastq} --extract
  """

}

process STAR{

  label 'big_mem'
  cpus 24
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_star/$animal_id/$hpi/$exp/$replicate"

  input:
  file indexes from star_indexes.collect()
  set animal_id, hpi, exp, replicate, complete_id,
      file(fastq) from fastqs_mapping

  output:
  set animal_id, hpi, exp, replicate, complete_id,
        file("${complete_id}.Log.final.out"),
        file("${complete_id}.Aligned.out.sam") into (mapped_sam, test)



  script:
  """
  STAR --genomeDir ${indexes} \
        --runThreadN ${task.cpus} \
        --readFilesIn ${fastq}  \
        --readFilesCommand zcat \
        --outFileNamePrefix ${complete_id}. \
        --outTmpDir $TMPDIR/${complete_id}.tmp
  """
}



process sort_bam{
  cpus 24
  label "rnaseq"
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_star/$animal_id/$hpi/$exp/$replicate"


  input:
  set animal_id, hpi, exp, replicate, complete_id, file(summary), file(sam) from mapped_sam

  output:
  set complete_id, animal_id, hpi, exp, replicate, \
      file(sam), file("${complete_id}_sorted.bam") into sorted_bams

  script:
  """
  picard-tools SortSam I=${sam} TMP_DIR=$TMPDIR/${complete_id}.tmp O=${complete_id}_sorted.bam SO=queryname
  """
}


sorted_bams.into{sorted_bams1; sorted_bams2; }

// process runRSeQC{
//
//   storeDir "${params.output_dir}/03_star/$animal_id/$hpi/$exp/$replicate/rseqc"
//
//   input:
//   set complete_id, animal_id, hpi, exp, replicate,
//       file(bam),
//       file(bai) from sorted_bams2
//   set file(bed), file(bed12) from  bed_channel.collect()
//
//   output:
//   file "*" into read_distribution_channel
//
//   script:
//   """
//   read_distribution.py -i ${bam} -r ${bed12} > ${complete_id}.UMI.f3.q60.read_distribution.txt
//   infer_experiment.py -i ${bam} -r ${bed12} > ${complete_id}.infer_experiment.txt
//   junction_annotation.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12}
//   bam_stat.py -i ${bam} > ${complete_id}.bam_stat.txt
//   junction_saturation.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12} > ${complete_id}.junction_annotation_log.txt
//   inner_distance.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12}
//   read_duplication.py -i ${bam} -o ${complete_id}.read_duplication
//   """
// }


process RevertSam{
  cpus 24
  label "rnaseq"
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_star/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id, file(bam) from bams_3

  output:
  set complete_id,file("${complete_id}_reverted.bam") into reverted_bams

  script:
  """
  picard-tools RevertSam I=${bam} O=${complete_id}_reverted.bam TMP_DIR=$TMPDIR/${complete_id}.tmp
  """
}


unmapped_and_mapped_bams = sorted_bams1.combine(reverted_bams, by:0)

unmapped_and_mapped_bams.into{unmapped_and_mapped_bams0; unmapped_and_mapped_bams1}

unmapped_and_mapped_bams1.subscribe{ println "$it" }

// // TODO: Missing orientation !
process MergeBamAlignment{

  cpus 12
  label "rnaseq"
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_star/$animal_id/$hpi/$exp/$replicate"

  input:
  file assembly from ReferenceChannel2.collect()
  file(dict) from DictChannel4.collect()
  set complete_id,animal_id, hpi, exp, replicate,
      file(sam), file(bam),
      file(unmapped_bam) from unmapped_and_mapped_bams0

  output:
  set animal_id, hpi, exp, replicate, complete_id, file("${complete_id}.merged.bam") into merged_bam_alignments_channel


  script:
  """
  picard-tools MergeBamAlignment REFERENCE_SEQUENCE=${assembly} \
               UNMAPPED_BAM=${unmapped_bam} \
               ALIGNED_BAM=${bam} \
               OUTPUT=${complete_id}.merged.bam \
               INCLUDE_SECONDARY_ALIGNMENTS=false \
               PAIRED_RUN=true \
               ORIENTATIONS=${params.strandness} \
               TMP_DIR=$TMPDIR/${complete_id}.tmp
  """
}

process TagReadWithGeneFunction{

  storeDir "${params.output_dir}/03_DropSeqPreProcessing/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id, file(bam) from merged_bam_alignments_channel
  file(gtf) from GtfChannel5.collect()

  output:
  set animal_id, hpi, exp, replicate, complete_id, file("${complete_id}_gene_exon_tagged.bam") into bams_tagged


  script:
  """
  TagReadWithGeneFunction I=${bam} \
                              O=${complete_id}_gene_exon_tagged.bam \
                              ANNOTATIONS_FILE=${gtf} \
                              USE_STRAND_INFO=true
  """

}
params.primer = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
process DetectBeadSubstitutionError{

  cpus 8
  storeDir "${params.output_dir}/03_DropSeqPreProcessing/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id, file(bam) from bams_tagged

  output:
  set animal_id, hpi, exp, replicate, complete_id, file("${complete_id}_clean.bam") into bams_tagged_clean
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

  cpus 8
  storeDir "${params.output_dir}/04_CellSelection/$animal_id/$hpi/$exp/$replicate"

  input:
  set animal_id, hpi, exp, replicate, complete_id, file(bam) from bams_tagged_clean

  output:
  set animal_id, hpi, exp, replicate, complete_id, file(bam), file("${complete_id}_out_cell_readcounts.txt.gz") into BamTagHistogramChannel


  script:
  """
  BamTagHistogram I=${bam} O=${complete_id}_out_cell_readcounts.txt.gz TAG=XC
  """
}

process CellSelection{

  cpus 4
  storeDir "${params.output_dir}/03_DropSeqPreProcessing/$animal_id/$hpi/$exp/$replicate"
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

// OUTPUT_READS_INSTEAD=Boolean  Output number of reads instead of number of unique molecular barcodes.  Default value:
// false. This option can be set to 'null' to clear the default value. Possible values:{true, false}
process DigitalExpressionMatrix_UMI{

  cpus 24
  storeDir "${params.output_dir}/04_DigitalExpressionMatrix/$animal_id/$hpi/$exp/$replicate"

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
