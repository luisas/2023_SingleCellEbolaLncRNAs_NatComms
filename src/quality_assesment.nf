#!/usr/bin/env nextflow

/*
 * Nextflow pipeline for quality assessment of RNA-seq data.
 */

// ------------  INPUT PARAMETERS -------------
params.dataset_bam_dir = "/gpfs/projects/bsc83/Ebola/00_RawData/pardis_shared_data/sabeti-txnomics/alin/190713_Zyagen-longRNA/tmp/00_demux/bams_per_lane/*/"
//params.dataset_bam_dir = "/gpfs/scratch/bsc83/bsc83024/test_dataset/bams_per_lane/*/"

// Folder where the output directories of the pipeline will be placed
params.output_dir = "/gpfs/projects/bsc83/Ebola/data/"
//params.output_dir = "/gpfs/scratch/bsc83/bsc83024/test_output_last/"

// If set to true the quality assessment will be computed with fastqc
params.fastqc = "true"
// params needed for mapping
params.assembly = "/gpfs/projects/bsc83/Ebola/00_InformationFiles/indexes/hisat2/rheMac8_EBOV-Kikwit"
params.assembly_prefix = "${params.assembly}".tokenize('/')[-1]
// Annotation
params.ss = "/gpfs/projects/bsc83/Ebola/00_InformationFiles/gene_annotations/rheMac8_EBOV-Kikwit.ss.txt"




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

// Channel for the Assembly
assemblyForMapping = Channel.fromPath("${params.assembly}*")
// Channel for Annotation
annotationForMapping = Channel.fromPath("${params.ss}")




/*
* STEP 1: Convert the unmapped bam files into fastq files.
* Raw reads are stored in the bam format but are not mapped yet. They need
* to be converted into fastq format to proceed with the pipeline.
*/
process convert_bam_to_fastq {

    storeDir "${params.output_dir}/01_fastq/$dataset_name/$tissue/$sample"

    input:
    set lane,lane_number, dataset_name, dayPostInfection, tissue, sample,complete_id, file(unmapped_bam) from dataset_bam

    output:
    set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
        file("${complete_id}.1.fq"),
        file("${complete_id}.2.fq"),
        file("${complete_id}.unpaired.fq") into (fastq_files_for_qc, fastq_files_for_mapping)

    script:
    """
    java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar SamToFastq I=${unmapped_bam} FASTQ=${complete_id}.1.fq  SECOND_END_FASTQ=${complete_id}.2.fq UNPAIRED_FASTQ=${complete_id}.unpaired.fq
    """

}

/*
* STEP 2.1: Generate quality assessment with fastqc.
*/
process generate_fastqc{

  storeDir "${params.output_dir}/02_fastqc/$dataset_name/$tissue/$sample"

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


/*
* STEP 2.2: Map, sort and index
*/
process mapping_hisat{

  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$sample"

  input:
  file assembly from assemblyForMapping.collect()
  file ss from annotationForMapping.collect()
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(fastq_1),file(fastq_2),file(unpaired)  from fastq_files_for_mapping

  output:
  set file("${complete_id}.novel_ss.txt"), file("${complete_id}.hisat2_summary.txt"),
      file("${complete_id}.sam"), file("${complete_id}.bam") into mapped_bams

  /*
  *  --known-splicesite-infile <path>   provide a list of known splice sites
  *  --novel-splicesite-outfile <path>  report a list of splice sites
  *  --dta                              reports alignments tailored for transcript assemblers
  *  --summary-file                     Print alignment summary to this file.
  */

  script:
  """
  hisat2-align-s -x ${params.assembly_prefix} -1  ${fastq_1} -2 ${fastq_2} --known-splicesite-infile ${ss} --novel-splicesite-outfile ${complete_id}.novel_ss.txt --downstream-transcriptome-assembly  --time --summary-file ${complete_id}.hisat2_summary.txt --rna-strandness FR > ${complete_id}.sam
  samtools sort -T $TMPDIR/${complete_id}.tmp  -O bam -o ${complete_id}.bam ${complete_id}.sam
  samtools index ${complete_id}.bam
  """

}
