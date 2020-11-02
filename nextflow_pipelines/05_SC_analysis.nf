// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"

params.output_dir_name = "01_scRNA-Seq_inVivo_rhemac10"
params.primer = "AAGCAGTGGTATCAACGCAGAGTGAATGGG"
params.strandness = "FR"
params.dataset_bam_dir = "${params.dirData}/00_RawData/seqwell/data/IRF_SerialSac/Mapping_V4/Mapping_V4"

params.output_dir_preliminary = "${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/"
params.output_dir = "${params.dirData}/02_scRNA-Seq_PBMCs/${params.output_dir_name}"


Channel.fromPath("${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")
       .into{ GtfChannel;GtfChannel2; GtfChannel3;GtfChannel4; GtfChannel5;   }

if("${params.output_dir_name}" == "01_scRNA-Seq_inVivo_rhemac10"){
  Channel.fromFilePairs("${params.dataset_bam_dir}/*/results/samples/*/trimmmed_repaired_R{1,2}.fastq.gz")
                .ifEmpty('bam files directory is empty')
                .map{tuple(it[1][0].toString().split('/')[-2].split('\\.')[0],
                           it[1][0].toString().split('/')[-2].split('\\.')[1],
                          it[1][0].toString().split('/')[-2].split('\\.')[2],
                          it[1][0].toString().split('/')[-2].split('\\.')[3],
                          it[1][0].toString().split('/')[-2].split('\\.')[4],
                          it[1][0].toString().split('/')[-2].split('\\.').join('_'),
                           it[1])}
                .into{ fastqs; fastqs2 }

}else if("${params.output_dir_name}" == "00_scRNA-Seq_exVivo_rhemac10"){

  Channel.fromFilePairs("${params.dataset_bam_dir}/20190829_Novaseq_88bp_ExVivo2*/*/E*_R{1,2}.fastq.gz")
          .ifEmpty('bam files directory is empty')
          .map{tuple(it[0].toString().split('\\.')[0],
                     it[0].toString().split('\\.')[1],
                    it[0].toString().split('\\.')[3],
                    it[0].toString().split('\\.')[6],
                    it[0].toString().split('\\.')[4],
                    it[0].toString().split('\\.')[0]+ "_" + it[0].toString().split('\\.')[6] + "_" + it[0].toString().split('\\.')[1] + "_" + it[0].toString().split('\\.')[3],
                     it[1])}
          .into{ fastqs; fastqs2; printing }
}
params.scripts="${baseDir}/scripts/"
gtfToGenePred_script_ch = Channel
                      .fromPath("${params.scripts}/gtfToGenePred")
genePredToBed_script_ch = Channel
                      .fromPath("${params.scripts}/genePredToBed")
cell_selection_script = Channel.fromPath("${params.scripts}/cellselection.R").collect()

Channel.fromPath("${params.output_dir_preliminary}/gene_annotations/rheMac10_EBOV-Kikwit.gtf")
        .into{ GtfNotNovel;   }

Channel.fromPath("${params.output_dir_preliminary}/reference_assembly/rheMac10_EBOV-Kikwit.fa")
       .into{ ReferenceChannel;ReferenceChannel2;   }

Channel.fromPath("${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/reference_assembly/rheMac10_EBOV-Kikwit.dict")
       .into{ DictChannel; DictChannel2;DictChannel3; DictChannel4;  }



log.info "=============================================="
log.info " Metadata Generation  "
log.info "=============================================="


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




process create_star_indexes{

  storeDir "${params.output_dir_preliminary}/indexes/"
  //cpus 16
  cpus 48
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

      process STAR{

        label 'big_mem'
        cpus 48
        //cpus 1
        tag "${complete_id}"
        storeDir "${params.output_dir}/02_star/$animal_id/$hpi/$exp/$replicate/$preprocessing"

        input:
        file indexes from star_indexes.collect()
        set animal_id, hpi, exp, replicate, preprocessing, complete_id,
            file(fastq_pair) from fastqs

        output:
        set animal_id, hpi, exp, replicate, preprocessing, complete_id,
              file("${complete_id}.Log.final.out"),
              file("${complete_id}.Aligned.out.bam") into (mapped_bam, test)



        script:
        """
        STAR --genomeDir ${indexes} \
              --runThreadN ${task.cpus} \
              --readFilesIn ${fastq_pair[0]} ${fastq_pair[1]}   \
              --readFilesCommand zcat \
              --outFileNamePrefix ${complete_id}. \
              --outTmpDir $TMPDIR/${complete_id} \
              --outSAMtype BAM Unsorted
        """
      }



process sort_bam{
  cpus 48
  //cpus 1
  label "rnaseq"
  tag "${complete_id}"
  storeDir "${params.output_dir}/02_star/$animal_id/$hpi/$exp/$replicate/$preprocessing"

  input:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file(summary), file(bam) from mapped_bam

  output:
  set complete_id, animal_id, hpi, exp, replicate, preprocessing, \
      file(bam), file("${complete_id}_sorted.bam") into sorted_bams

  script:
  """
  picard-tools SortSam I=${bam} O=${complete_id}_sorted.bam SO=queryname TMP_DIR=$PWD/tmp
  """
}


sorted_bams.into{sorted_bams1; sorted_bams2; }

process runRSeQC{

  storeDir "${params.output_dir}/02_star/$animal_id/$hpi/$exp/$replicate/$preprocessing"
  input:
  set complete_id, animal_id, hpi, exp, replicate, preprocessing,
    file(bam),
    file(sortedbam) from sorted_bams2
    set file(bed), file(bed12) from  bed_channel.collect()

  output:
   file "*" into read_distribution_channel

   script:
   """
   infer_experiment.py -i ${sortedbam} -r ${bed12} > ${complete_id}.infer_experiment.txt
   """
}


process FastqToSam{
  cpus 48
  //cpus 1
  label "rnaseq"
  tag "${complete_id}"
  storeDir "${params.output_dir}/02_star/$animal_id/$hpi/$exp/$replicate/$preprocessing"
  input:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id,
      file(fastq_pair) from fastqs2
  output:
  set complete_id,file("${complete_id}_unmapped.bam") into reverted_bams

  script:
  """
  picard-tools FastqToSam FASTQ=${fastq_pair[0]} FASTQ2=${fastq_pair[1]} O=${complete_id}_unmapped.bam TMP_DIR=$TMPDIR/${complete_id}.tmp SAMPLE_NAME=${complete_id}
  """
}


unmapped_and_mapped_bams = sorted_bams1.combine(reverted_bams, by:0)

unmapped_and_mapped_bams.into{unmapped_and_mapped_bams0; unmapped_and_mapped_bams1}

unmapped_and_mapped_bams1.subscribe{ println "$it" }

// // TODO: Missing orientation !
process MergeBamAlignment{

  cpus 48
  //cpus 1
  label "rnaseq"
  tag "${complete_id}"
  storeDir "${params.output_dir}/02_star/$animal_id/$hpi/$exp/$replicate/$preprocessing"
  input:
  file assembly from ReferenceChannel2.collect()
  file(dict) from DictChannel4.collect()
  set complete_id,animal_id, hpi, exp, replicate, preprocessing,
      file(sam), file(bam),
      file(unmapped_bam) from unmapped_and_mapped_bams0

  output:
  set animal_id, hpi, exp, replicate, preprocessing, complete_id, file("${complete_id}.merged.bam") into merged_bam_alignments_channel


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
