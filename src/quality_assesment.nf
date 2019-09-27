#!/usr/bin/env nextflow

/*
 * TODO COPYRIGHT
 */

 log.info "=============================================="
 log.info "RNASeq analysis pipeline for Zyagen Samples"
 log.info "=============================================="

 /* --------------------------------------
  *  INPUT PARAMETERS defaults
  * --------------------------------------
  */
// Unmapped bam files folder - RAW DATA
params.dataset_bam_dir = "/gpfs/projects/bsc83/Data/Ebola/00_RawData/pardis_shared_data/sabeti-txnomics/alin/190713_Zyagen-longRNA/tmp/00_demux/bams_per_lane/*/"
//params.dataset_bam_dir = "/gpfs/scratch/bsc83/bsc83024/test_dataset/bams_per_lane/*/"

// Folder where the output directories of the pipeline will be placed
params.output_dir = "/gpfs/projects/bsc83/Ebola/data/02_RNA-Seq/"
//params.output_dir = "/gpfs/projects/bsc83/Ebola/test_data/02_test/"


params.preliminary_files_dir="/gpfs/projects/bsc83/Ebola/data/01_PreliminaryFiles"
params.assembly_name = "rheMac8_EBOV-Kikwit"

// This ones should be automatically retrieved but can be changed in need
params.reference_assembly = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.fa"
// Dictionary
params.dict = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.dict"
// Indexes
params.hisat2_indexes = "${params.preliminary_files_dir}/indexes/hisat2/${params.assembly_name}"
// Known splice sites
params.ss = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.ss.txt"
// Gene annotation
params.gene_annotation = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.bed"

// If set to true the quality assessment will be computed with fastqc
params.fastqc=true

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
* - complete_id (${dataset_name}_${dayPostInfection}_${tissue}_${sample}_l${lane_number})
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
                           it.baseName.split('_')[0] + "_" + it.baseName.split('_')[2].split('-')[0] +"_"+ it.baseName.split('_')[1] +"_"+ it.baseName.split('_')[2].split('-')[1] +"_"+ "l" +  it.parent.name.split('\\.')[1],
                           it ) }


// the unmapped_bam channel contains the exact same files as the dataset_bam one.
// It only differs in the values that it maps to it, in this case only the complete_id.
// It is used for the add_unmapped_bam process.
unmapped_bams = Channel
                .fromPath("${params.dataset_bam_dir}/*_long.bam")
                .ifEmpty('bam files directory is empty')
                .map{ tuple(it.baseName.split('_')[0] + "_" + it.baseName.split('_')[1] +"_"+ it.baseName.split('_')[2].split('-')[0] +"_"+ it.baseName.split('_')[2].split('-')[1] +"_"+ "l" +  it.parent.name.split('\\.')[1],
                            it)}

// Channel for the Indexes
indexesForMapping = Channel.fromPath("${params.hisat2_indexes}*")
//Channel for the reference assembly
reference_assembly_channel = Channel.fromPath("${params.reference_assembly}*")
// Channel for known splice sites
known_ss = Channel.fromPath("${params.ss}")

Channel.fromPath("${params.gene_annotation}")
                                 .into{ gene_annotation_channel; gene_annotation_channel_2}
dictionary_channel = Channel.fromPath("${params.dict}")


/*  ----------------------------------------------------------------------
*   ----------------------------------------------------------------------
*                       BEGINNING OF THE PIPELINE
*   ----------------------------------------------------------------------
*   ----------------------------------------------------------------------
*/


/*
* STEP 1: Convert the unmapped bam files into fastq files.
* Raw reads are stored in the bam format but are not mapped yet. They need
* to be converted into fastq format to proceed with the pipeline.
*/
process convert_bam_to_fastq {

    tag "${complete_id}"
    //publishDir  "${params.output_dir}/01_fastq/$dataset_name/$tissue/$dayPostInfection/$sample"
    publishDir "${params.output_dir}/01_fastq/$dataset_name/$tissue/$dayPostInfection/$sample" , mode: 'copy'

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
* STEP 2: Generate quality assessment with fastqc.
*/
process generate_fastqc{

  tag "${complete_id}"
  publishDir "${params.output_dir}/02_fastqc/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'

  when:
  params.fastqc

  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(fastq_1),file(fastq_2),file(unpaired)  from fastq_files_for_qc

  output:
  file "*" into fastqcs

  script:
  """
  mkdir -p output_fastq_dir
  fastqc ${fastq_1} ${fastq_2} --extract
  """

}


/*
* STEP 3.1: Mapping with HISAT.
* It maps the reads of the generated fastq files
* to the assembly defined in params.assembly
* The indexes were computed WITHOUT the help of known splice sites
* the known splice sites are now used as input for the mapper.
*/
process mapping_hisat{

  label 'big_mem'
  tag "${complete_id}"
  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample" , mode: 'copy'

  input:
  file indexes from indexesForMapping.collect()
  file ss from known_ss.collect()
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(fastq_1),file(fastq_2),file(unpaired)  from fastq_files_for_mapping

  output:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file("${complete_id}.novel_ss.txt"), file("${complete_id}.hisat2_summary.txt"),
      file("${complete_id}.sam") into mapped_sam

  /*
  *  --known-splicesite-infile <path>   provide a list of known splice sites
  *  --novel-splicesite-outfile <path>  report a list of splice sites
  *  --dta                              reports alignments tailored for transcript assemblers
  *  --summary-file                     Print alignment summary to this file.
  */

  script:
  """
  hisat2-align-s -x ${params.assembly_name} -1  ${fastq_1} -2 ${fastq_2} --known-splicesite-infile ${ss} --novel-splicesite-outfile ${complete_id}.novel_ss.txt --downstream-transcriptome-assembly  --time --summary-file ${complete_id}.hisat2_summary.txt --rna-strandness FR > ${complete_id}.sam
  """

}

/*
* STEP 3.2: Sort and index the mapped bam files
*
* The sorting is done with picard as we encountered some problems when running
* samtools and afterwards other picardTools on top of the sorted bam.
* TODO add the indexing
*/
process sort_and_index_bam{

  tag "${complete_id}"
  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'


  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(ss), file(summary),
      file(sam) from mapped_sam

  output:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(ss), file(summary),
      file(sam), file("${complete_id}.bam"), file("${complete_id}.bam.bai") into (sorted_and_index_bam, sorted_indexed_bams_for_stats)


  script:
  """
  samtools sort ${sam} -T ${complete_id}.tmp  -O bam -o ${complete_id}.bam
  samtools index ${complete_id}.bam
  """
}

//java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar SortSam INPUT=${sam} OUTPUT=${complete_id}.bam SORT_ORDER=cordinate
//java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar BuildBamIndex I=${complete_id}.bam

/*
* STEP 4: Merge bam alignments -  Add UMIs
*
* After having converted the unmapped bam files into fastq we have lost
* relevant information like UMIs info.
* We now merge the unmapped and the mapped bam files to retrieve information
* about the UMIs.
*/

// We join the mapped and unmapped channel by complete_id.
// See how channel were built at the beginning - complete id is always in the -2 position
unmapped_and_mapped_bams = sorted_and_index_bam.join(unmapped_bams, by:-2)

process add_unmapped_bam{

  tag "${complete_id}"
  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'

  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(ss), file(summary),
      file(sam), file(bam),file(bai),
      complete_id_unmapped, file(unmapped_bam) from unmapped_and_mapped_bams
  file assembly from reference_assembly_channel.collect()
  file dict from dictionary_channel.collect()

  output:
  file ("${complete_id}.UMI.bam") into merged_bam_alignments_channel

  script:
  """
  # sort and index the unmapped_bam
  samtools sort ${unmapped_bam} -T ${complete_id}_long.tmp  -O bam -o ${complete_id}_unmapped.bam
  samtools index ${complete_id}_unmapped.bam
  java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar MergeBamAlignment UNMAPPED=${complete_id}_unmapped.bam ALIGNED=${bam} O=${complete_id}.UMI.bam R=${assembly} SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
  """
}

//java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar SortSam INPUT=${unmapped_bam} OUTPUT=${complete_id_unmapped}-long-sorted.bam SORT_ORDER=coordinate
//java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar BuildBamIndex -I=${complete_id_unmapped}-long-sorted.bam




// In this channel we group the files by the complete_id after removing the lane.
// We have now a set of bams per sample from all the lanes.
merged_bam_alignments_channel.map { file ->
        def key = file.name.toString().tokenize('_').init()
        key_string = key.join("_")
        key << key_string
        return tuple(key, file)
     }
     .groupTuple()
     .map{ tuple(it.get(0).get(0),
                 it.get(0).get(1),
                 it.get(0).get(2),
                 it.get(0).get(3),
                 it.get(0).get(4),
                 it.get(1))}
    .set{ groups_lanes }

/*
* Step 5 Merge lanes
* Each sample has 4 lanes.
* We need to merge all of them into one.
*
* Aaron's mail
* We prepared both the old and new Zyagen samples at the same time with the s
* same library prep. We introduce UMIs at the adaptor ligation step, prior to any PCR.
* The bam files in the recent run folder represent a single library prep per
* sample, sequenced over four lanes, for both the old and new Zyagen samples.
* So we should merge bams for all four lanes, and then
*/
process merge_lanes{

  tag "${complete_id}"
  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(samples) from groups_lanes

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file("${complete_id}.UMI.bam"), file("${complete_id}.UMI.bam.bai") into merged_bylanes

  script:
  """
  samtools merge -f ${complete_id}.UMI.bam  ${samples}
  # Merge mapped bam files with UMI info
  samtools index  ${complete_id}.UMI.bam
  """
}


// Filter put bams - now it is stringent

/*
*   STEP 6 filter bams
*   Filter low-quality reads and retain properly-paired uniquely mapped reads
*   f3 : read paired & read mapped in proper pair
*   q60: 60 should be mapped TODO check properly
*   Now it is a stringent filtering. We may wanna make it looser for the lncRNA
*   annotation pipeline.
*/
process filter_bams_samtools{

  tag "${complete_id}"
  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(merged_bam), file(bai_mereged_bam) from merged_bylanes

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file("${complete_id}.UMI.f3.q60.bam"), file("${complete_id}.UMI.f3.q60.bam.bai") into (filtered_merged_bams_1, filtered_merged_bams_2)

  script:
  """
  samtools view -b -f3 -q 60  ${merged_bam} > ${complete_id}.UMI.f3.q60.bam
  samtools index ${complete_id}.UMI.f3.q60.bam
  """
}

/*
* STEP 7A : Mark duplicates with picard
* we probably need to delete it.
*/

process MarkDuplicates_picard{

  tag "${complete_id}"
  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(filtered_bam), file(filtered_bai) from filtered_merged_bams_1

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,file("${complete_id}.md.bam"), file("${complete_id}.md_metrics.bam") into marked_d_bams

  script:
  """
  java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar MarkDuplicates I=${filtered_bam} O=${complete_id}.UMI.f3.q60.md.bam M=${complete_id}.UMI.f3.q60.md_metrics.bam REMOVE_DUPLICATES=false ASSUME_SORTED=true  CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
  """
}


/*
* STEP 7B : Dedup duplicate using UMIs
*/

process dedupUmi{

  tag "${complete_id}"
  label 'big_mem'
  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(filtered_bam), file(filtered_bai) from filtered_merged_bams_2

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file("${complete_id}.UMI.f3.q60.umi_dedup.bam"),
      file("${complete_id}.UMI.f3.q60.umi_dedup.bam.bai") into (dedup_umi_bams, dedup_umi_bams_for_count, dedup_umi_bams_for_distr)

  script:
  """
  umi_tools dedup --stdin=${filtered_bam} --log=${complete_id}.umi_tools.dedup.log --output-stats=${complete_id} --extract-umi-method=tag --umi-tag=RX --paired --method directional  > ${complete_id}.UMI.f3.q60.umi_dedup.bam
  samtools index ${complete_id}.UMI.f3.q60.umi_dedup.bam
  """

}

/*
*  STEP 8.1 count reads per chromosome to grasp an overview.
*/
process getCountsUMIs{

  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(umi_bam),
      file(umi_bai) from dedup_umi_bams_for_count

  output:
  file "counts/*" into counts_umis

  script:
  """
  mkdir counts
  # Create all the chrom names
  chroms=""
  for i in {1..20};do chroms=\${chroms}" chr"\${i};done

  n=`samtools view ${umi_bam} -c`
  n_autosomes=`samtools view  ${umi_bam} \${chroms} -c`
  n_chrX=`samtools view  ${umi_bam} chrX -c`
  n_chrY=`samtools view  ${umi_bam} chrY -c`
  n_chrM=`samtools view  ${umi_bam} chrM -c`
  n_EBOV=`samtools view  ${umi_bam} EBOV_Kikwit -c`
  n_chrUn=`samtools view ${umi_bam} | awk '{if(\$3~/chrUn/)print}' | wc -l | awk '{print \$1}'`

  echo -e "Total\t"\${n} > counts/${complete_id}.UMI.f3.q60.n_reads.txt
  echo -e "Autosomes\t"\${n_autosomes} >> counts/${complete_id}.UMI.f3.q60.n_reads.txt
  echo -e "chrX\t"\${n_chrX} >> counts/${complete_id}.UMI.f3.q60.n_reads.txt
  echo -e "chrY\t"\${n_chrY} >> counts/${complete_id}.UMI.f3.q60.n_reads.txt
  echo -e "chrM\t"\${n_chrM} >> counts/${complete_id}.UMI.f3.q60.n_reads.txt
  echo -e "EBOV\t"\${n_EBOV} >> counts/${complete_id}.UMI.f3.q60.n_reads.txt
  """
}


/*
*  STEP 8.2 get read distribution numbers
*/

process getReadDistribution{

  publishDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample", mode: 'copy'

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(umi_bam),
      file(umi_bai) from dedup_umi_bams_for_distr
  file(gene_annotation) from gene_annotation_channel_2.collect()

  output:
  file "read_distribution/*" into read_distribution_channel

  script:
  """
  mkdir read_distribution
  read_distribution.py -i ${umi_bam} -r ${gene_annotation} > read_distribution/${complete_id}.UMI.f3.q60.read_distribution.txt
  """
}
// -----------------------------------------------------------------------------


/*   -------------------------------
*           Groovy Functions
*    -------------------------------
*/

def remove_lane_from_id(String id){
  return id.split("_").init().join("_")
}
