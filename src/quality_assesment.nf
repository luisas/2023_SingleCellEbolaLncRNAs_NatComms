#!/usr/bin/env nextflow

/*
 * Nextflow pipeline for quality assessment of RNA-seq data.
 */

// ------------  INPUT PARAMETERS -------------
params.dataset_bam_dir = "/gpfs/scratch/bsc83/bsc83024/test_dataset/bams_per_lane/*/"

// Folder where the output directories of the pipeline will be placed
params.output_dir = "/gpfs/scratch/bsc83/bsc83024/test_output/"

// If set to true the quality assessment will be computed with fastqc
params.fastqc = "true"

// -----------------------------------------------

/* Channel for all the unmapped bam files present in any of the dataset_bam_dir
* subdirectories.
* In order, it maps:
* - lane
* - lane_number
* - dataset_name
* - dayPostInfection
* - tissue
* - sample
* - complete_id (${tissue}.${sample}.l${lane_number})
* - the file itself.
*/
dataset_bam = Channel
              .fromPath("${params.dataset_bam_dir}/*_long.bam")
              .ifEmpty('bam files directory is empty')
              .map { tuple(it.parent.name,
                           it.parent.name.split('\\.')[1],
                           it.baseName.split('_')[0],
                           it.baseName.split('_')[1],
                           it.baseName.split('_')[2].split('-')[0],
                           it.baseName.split('_')[2].split('-')[1],
                           it.baseName.split('_')[2].split('-')[0] + it.baseName.split('_')[2].split('-')[1] + "l." +  it.parent.name.split('\\.')[1],
                           it ) }

/*
* STEP 1: Convert the unmapped bam files into fastq files.
*/
process convert_bam_to_fastq {

    publishDir "${params.output_dir}/01_fastq/$dataset_name/$tissue/$sample"

    input:
    set lane,lane_number, dataset_name, dayPostInfection, tissue, sample,complete_id, file(unmapped_bam) from dataset_bam

    output:
    set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
        file("${complete_id}.1.fq"),
        file("${complete_id}.2.fq"),
        file("${complete_id}.unpaired.fq") into (fastq_files_for_qc, fastq_files)

    script:
    """
    java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar SamToFastq I=${unmapped_bam} FASTQ=${complete_id}.1.fq  SECOND_END_FASTQ=${complete_id}.2.fq UNPAIRED_FASTQ=${complete_id}.unpaired.fq
    """

}

/*
* STEP 2: Generate quality assessment with fastqc.
*/
process generate_fastqc{

  publishDir "${params.output_dir}/02_fastqc/$dataset_name/$tissue/$sample"

  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(fastq_1),file(fastq_2),file(unpaired)  from fastq_files_for_qc

  output:
  file "*" into fastqcs

  when:
  params.fastqc == "true"

  script:
  """
  mkdir -p output_fastq_dir
  fastqc ${fastq_1} ${fastq_2} --extract
  """

}
