#!/usr/bin/env nextflow

/*
 * TODO COPYRIGHT
 */

 log.info "=================================================================================="
 log.info "     RNASeq analysis pipeline & de Novo Transcriptome assembly"
 log.info "=================================================================================="

 /* --------------------------------------
  *  INPUT PARAMETERS defaults
  * --------------------------------------
  */
// Unmapped bam files folder - RAW DATA
params.dataset_bam_dir = "/gpfs/projects/bsc83/Data/Ebola/00_RawData/pardis_shared_data/sabeti-txnomics/alin/190713_Zyagen-longRNA/tmp/00_demux/bams_per_lane/*/"
//params.dataset_bam_dir = "/gpfs/scratch/bsc83/bsc83024/test_dataset/bams_per_lane/*/"
//params.dataset_bam_dir = "/gpfs/projects/bsc83/Data/Ebola/00_RawData_links/"

// Folder where the output directories of the pipeline will be placed
params.output_dir = "/gpfs/projects/bsc83/Projects/Ebola/data/02_RNA-Seq/"
//params.output_dir = "/gpfs/projects/bsc83/Projects/Ebola/test_data/04_test/"


params.preliminary_files_dir="/gpfs/projects/bsc83/Projects/Ebola/data/01_PreliminaryFiles"
params.assembly_name = "rheMac8_EBOV-Kikwit"

// This ones should be automatically retrieved but can be changed in need
params.reference_assembly = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.fa"
params.reference_assembly_fai = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.fa.fai"
// Dictionary
params.dict = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.dict"
// Indexes
params.hisat2_indexes = "${params.preliminary_files_dir}/indexes/hisat2/${params.assembly_name}"
// Known splice sites
params.ss = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.ss.txt"
// Gene annotation
params.gene_annotation = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.bed"
params.gtf = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.gtf"
params.gtf_rheMac = "${params.preliminary_files_dir}/gene_annotations/rheMac8.ensembl_release97.gtf"
params.gtf_transcriptonly = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}_rnaseqc.gtf"

params.bed_rheMac = "/gpfs/projects/bsc83/Data/Ebola/00_RawData/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8/Ensembl/rheMac8.Ensembl.bed"
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

params.zyagen = "false"
if("${params.zyagen}" != "false"){
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
                  .map{ tuple(it.baseName.split('_')[0] + "_" + it.baseName.split('_')[2].split('-')[0] +"_"+ it.baseName.split('_')[1] +"_"+ it.baseName.split('_')[2].split('-')[1] +"_"+ "l" +  it.parent.name.split('\\.')[1],
                              it)}
}else{
    dataset_bam = Channel
                  .fromPath("${params.dataset_bam_dir}*.bam", followLinks:true)
                  .map { tuple(it.baseName.split('_')[4].split('\\.')[0],
                               it.baseName.split('_')[4].split('\\.')[1],
                               it.baseName.split('_')[0],
                               it.baseName.split('_')[2],
                               it.baseName.split('_')[1],
                               it.baseName.split('_')[3],
                               it.baseName.split('_')[0] + "_" + it.baseName.split('_')[1]+ "_" +it.baseName.split('_')[2] + "_" +it.baseName.split('_')[3] + "_" + "l"+it.baseName.split('_')[4].split('\\.')[1],
                               it ) }
    // the unmapped_bam channel contains the exact same files as the dataset_bam one.
    // It only differs in the values that it maps to it, in this case only the complete_id.
    // It is used for the add_unmapped_bam process.
    unmapped_bams = Channel
                    .fromPath("${params.dataset_bam_dir}/*.bam")
                    .ifEmpty('bam files directory is empty')
                    .map{ tuple(it.baseName.split('_')[0] + "_" + it.baseName.split('_')[1]+ "_" +it.baseName.split('_')[2] + "_" +it.baseName.split('_')[3] + "_" + "l"+it.baseName.split('_')[4].split('\\.')[1],
                                it)}
}




// Channel for the Indexes
indexesForMapping = Channel.fromPath("${params.hisat2_indexes}*")
indexesForMapping2 = Channel.fromPath("${params.hisat2_indexes}*")

//Channel for the reference assembly
reference_assembly_channel = Channel.fromPath("${params.reference_assembly}*")
Channel.fromPath("${params.reference_assembly}").into{reference_genome_ch;reference_genome_ch_2; reference_genome_ch_3}
reference_genome_fai_ch = Channel.fromPath("${params.reference_assembly_fai}")
// Channel for known splice sites
known_ss = Channel.fromPath("${params.ss}")
known_ss_2 = Channel.fromPath("${params.ss}")
Channel.fromPath("${params.gene_annotation}")
                                 .into{ gene_annotation_channel; gene_annotation_channel_2; gene_annotation_channel_3}
dictionary_channel = Channel.fromPath("${params.dict}")
Channel.fromPath("${params.gtf_transcriptonly}").into{ gtf_transcriptonly_channel; gtf_transcriptonly_channel_2; gtf_transcriptonly_channel_3}
Channel.fromPath("${params.gtf}").into{ gtfChannel; gtfChannel2; gtfChannel3; gtfChannel4; gtfChannel5}
Channel.fromPath("${params.bed_rheMac}").into{ gene_annotation_rheMac_bed; gene_annotation_rheMac_bed_2}
Channel.fromPath("${params.gtf_rheMac}").into{ gtf_rheMac_channel; gtf_rheMac_channel2; gtf_rheMac_channel3}
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
    //storeDir  "${params.output_dir}/01_fastq/$dataset_name/$tissue/$dayPostInfection/$sample"
    storeDir "${params.output_dir}/01_fastq/$dataset_name/$tissue/$dayPostInfection/$sample"

    input:
    set lane,lane_number, dataset_name, dayPostInfection, tissue, sample,complete_id, file(unmapped_bam) from dataset_bam

    output:
    set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
        file("${complete_id}.1.fq"),
        file("${complete_id}.2.fq"),
        file("${complete_id}.unpaired.fq") into (fastq_files_for_qc, fastq_files_for_mapping, fastqfiles_for_trimmomatic)

    script:
    """
    java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar SamToFastq I=${unmapped_bam} FASTQ=${complete_id}.1.fq  SECOND_END_FASTQ=${complete_id}.2.fq UNPAIRED_FASTQ=${complete_id}.unpaired.fq
    """

}

// Remove adapters (ILLUMINACLIP:TruSeq3-PE.fa:2:30:10)
// Remove leading low quality or N bases (below quality 3) (LEADING:3)
// Remove trailing low quality or N bases (below quality 3) (TRAILING:3)
// Scan the read with a 4-base wide sliding window, cutting when the average quality per base drops below 15 (SLIDINGWINDOW:4:15)
// Drop reads below the 36 bases long (MINLEN:36)

process trimmomatic{

  tag "${complete_id}"
  storeDir "${params.output_dir}/01b_qualityFiltering/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(readOne),
      file(readTwo),
      file(unpaired) from fastqfiles_for_trimmomatic

  output:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file("${complete_id}.1.fq.gz"),
      file("${complete_id}.2.fq.gz"),
      file("${complete_id}_forward_unpaired.fq.gz") into (fastqfiles_trimmomatic_for_qc,fastqfiles_trimmomatic)

  script:
  """
  java -jar /apps/TRIMMOMATIC/0.39/trimmomatic-0.39.jar PE -phred33 ${readOne} ${readTwo} ${complete_id}.1.fq.gz ${complete_id}_forward_unpaired.fq.gz ${complete_id}.2.fq.gz ${complete_id}_reverse_unpaired.fq.gz LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36
  """
}

/*
* STEP 2: Generate quality assessment with fastqc.
*/
process generate_fastqc{

  tag "${complete_id}"
  storeDir "${params.output_dir}/02_fastqc/$dataset_name/$tissue/$dayPostInfection/$sample"

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

process generate_fastqc_trimmomatic{

  tag "${complete_id}"
  storeDir "${params.output_dir}/02_fastqc_filtered/$dataset_name/$tissue/$dayPostInfection/$sample"

  when:
  params.fastqc

  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(fastq_1),file(fastq_2),file(unpaired)  from fastqfiles_trimmomatic_for_qc

  output:
  file "*" into fastqcs_timmomatic

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
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

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
* STEP 3.1: Mapping with HISAT.
* It maps the reads of the generated fastq files
* to the assembly defined in params.assembly
* The indexes were computed WITHOUT the help of known splice sites
* the known splice sites are now used as input for the mapper.
*/
process mapping_hisat_trimmomatic{

  cpus 24
  label 'big_mem'
  tag "${complete_id}"
  storeDir "${params.output_dir}/03b_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  file indexes from indexesForMapping2.collect()
  file ss from known_ss_2.collect()
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(fastq_1),file(fastq_2),file(unpaired)  from fastqfiles_trimmomatic

  output:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file("${complete_id}.novel_ss.txt"), file("${complete_id}.hisat2_summary.txt"),
      file("${complete_id}.sam") into mapped_sam_trimmomatic

  /*
  *  --known-splicesite-infile <path>   provide a list of known splice sites
  *  --novel-splicesite-outfile <path>  report a list of splice sites
  *  --dta                              reports alignments tailored for transcript assemblers
  *  --summary-file                     Print alignment summary to this file.
  */

  script:
  """
  /gpfs/projects/bsc83/utils/hisat2-2.1.0/hisat2 -p ${task.cpus} -x ${params.assembly_name} -1  ${fastq_1} -2 ${fastq_2} --known-splicesite-infile ${ss} --novel-splicesite-outfile ${complete_id}.novel_ss.txt --downstream-transcriptome-assembly  --time --summary-file ${complete_id}.hisat2_summary.txt --rna-strandness FR > ${complete_id}.sam
  """

}
/*
* STEP 3.2: Sort and index the mapped bam files
*
* The sorting is done with picard as we encountered some problems when running
* samtools and afterwards other picardTools on top of the sorted bam.
*/
process sort_bam{

  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"


  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(ss), file(summary),
      file(sam) from mapped_sam

  output:
  set complete_id,lane,lane_number, dataset_name, dayPostInfection, tissue, sample,
      file(ss), file(summary),
      file(sam), file("${complete_id}.bam")   into (sorted_and_index_bam, sorted_indexed_bams_for_stats, hisat2_bams)

  script:
  """
  samtools sort ${sam} -T $TMPDIR/${complete_id}.tmp  -O bam -o ${complete_id}.bam
  """
}

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
unmapped_and_mapped_bams = sorted_and_index_bam.combine(unmapped_bams, by:0)

process add_unmapped_bam{

  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set complete_id,lane,lane_number, dataset_name, dayPostInfection, tissue, sample,
      file(ss), file(summary),
      file(sam), file(bam),
      file(unmapped_bam) from unmapped_and_mapped_bams
  file assembly from reference_assembly_channel.collect()
  file dict from dictionary_channel.collect()

  output:
  file ("${complete_id}.UMI.bam") into merged_bam_alignments_channel

  script:
  """
  java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar MergeBamAlignment UNMAPPED=${unmapped_bam} ALIGNED=${bam} O=${complete_id}.UMI.bam R=${assembly} SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
  """
}

// // In this channel we group the files by the complete_id after removing the lane.
// // We have now a set of bams per sample from all the lanes.
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
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, tissue, dayPostInfection, sample, complete_id, file(samples) from groups_lanes

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file("${complete_id}.UMI.bam"), file("${complete_id}.UMI.bam.bai") into (merged_bylanes, raw_bams)

  script:
  """
  samtools merge -f ${complete_id}.UMI.bam  ${samples}
  # Merge mapped bam files with UMI info
  samtools index  ${complete_id}.UMI.bam
  """
}



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
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(merged_bam), file(bai_mereged_bam) from merged_bylanes

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file("${complete_id}.UMI.f3.q60.bam"), file("${complete_id}.UMI.f3.q60.bam.bai") into (filtered_merged_bams_1, filtered_merged_bams_2,filtered_merged_bams_3, filtered_merged_bams_4, filtered_merged_bams_5)

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

process removeDuplicates_picard{

  cpus 4
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(filtered_bam), file(filtered_bai) from filtered_merged_bams_1

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,file("${complete_id}.UMI.f3.q60.md.bam"), file("${complete_id}.UMI.f3.q60.md.bam.bai") into  (marked_d_bams, md_bams, md_bams_3)

  script:
  """
  java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar MarkDuplicates I=${filtered_bam} O=${complete_id}.UMI.f3.q60.md.bam M=${complete_id}.UMI.f3.q60.md_metrics.bam REMOVE_DUPLICATES=true ASSUME_SORTED=true  CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
  samtools index ${complete_id}.UMI.f3.q60.md.bam
  """
}


/*
* STEP 7B : Dedup duplicate using UMIs
*/

process dedupUmi{

  tag "${complete_id}"
  label 'big_mem'
  cpus 8
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(filtered_bam), file(filtered_bai) from filtered_merged_bams_2

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file("${complete_id}.UMI.f3.q60.umi_dedup.bam"),
      file("${complete_id}.UMI.f3.q60.umi_dedup.bam.bai") into (dedup_umi_bams, dedup_umi_bams_for_count, dedup_umi_bams_for_distr, dedup_umi_bams_for_count_2, dedup_umis_3, dedup_umis_4)

  script:
  """
  umi_tools dedup --stdin=${filtered_bam} --log=${complete_id}.umi_tools.dedup.log --output-stats=${complete_id} --extract-umi-method=tag --umi-tag=RX --paired --method directional  > ${complete_id}.UMI.f3.q60.umi_dedup.bam
  samtools index ${complete_id}.UMI.f3.q60.umi_dedup.bam
  """

}

/*
*  Compute the HTseq Counts for all the
*/
process getHTseqCountsMD{

  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,file(md_bam), file(md_metrics_bam) from marked_d_bams
  file gtf from gtfChannel.collect()

  output:
  file "${md_file_prefix}.HTseq.gene_counts.tab" into htseqCountsMd_channel

  script:
  md_file_prefix = get_file_name_no_extension(md_bam.name)
  """
  # Calc the counts for the md duplicates
  htseq-count -f bam -r pos -s yes -t gene -i gene_id ${md_bam} ${gtf} > ${md_file_prefix}.HTseq.gene_counts.tab
  """
}

process getHTseqCountsFilter{

  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(filtered_bam), file(filtered_bai) from filtered_merged_bams_3
  file gtf from gtfChannel3.collect()

  output:
  file "${filtered_file_prefix}.HTseq.gene_counts.tab" into htseqCountsFiltered_channel

  script:
  filtered_file_prefix = get_file_name_no_extension(filtered_bam.name)
  """
  # Calc the counts for the md duplicates
  htseq-count -f bam -r pos -s yes -t gene -i gene_id ${filtered_bam} ${gtf} > ${filtered_file_prefix}.HTseq.gene_counts.tab
  """
}

process getHTseqCountsUMI{

  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(umi_dedup_bam), file(umi_dedup_bai) from dedup_umi_bams_for_count
  file gtf from gtfChannel4.collect()

  output:
  file "${umi_dedup_file_prefix}.HTseq.gene_counts.tab" into htseqCountsUmi_channel

  script:
  umi_dedup_file_prefix = get_file_name_no_extension(umi_dedup_bam.name)
  """
  # Calc the counts for the umi_dedup
  htseq-count -f bam -r pos -s yes -t gene -i gene_id ${umi_dedup_bam} ${gtf} > ${umi_dedup_file_prefix}.HTseq.gene_counts.tab
  """
}


/*
*  STEP 8.1 count reads per chromosome to grasp an overview.
*/
process getCountsUMIs{

  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/counts"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bam),
      file(bai) from dedup_umi_bams_for_count_2

  output:
  file "*.n_reads.txt" into counts_umis

  script:
  file_prefix = get_file_name_no_extension(bam.name)
  template 'getcounts'
}

process getCountsFiltered{

  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/counts"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bam),
      file(bai) from filtered_merged_bams_4

  output:
  file "*.n_reads.txt" into counts_filtered

  script:
  file_prefix = get_file_name_no_extension(bam.name)
  template 'getcounts'
}

process getCountsMDs{

  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/counts"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bam),
      file(bai) from md_bams

  output:
  file "*.n_reads.txt" into counts_md

  script:
  file_prefix = get_file_name_no_extension(bam.name)
  template 'getcounts'

}

/*
*  STEP 8.2 get read distribution numbers
*/

process getReadDistribution{

  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/read_distribution"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(umi_bam),
      file(umi_bai) from dedup_umi_bams_for_distr
  file(gene_annotation) from gene_annotation_rheMac_bed.collect()

  output:
  file "${complete_id}.UMI.f3.q60.read_distribution.txt" into read_distribution_channel

  script:
  """
  read_distribution.py -i ${umi_bam} -r ${gene_annotation} > ${complete_id}.UMI.f3.q60.read_distribution.txt
  """
}


raw_and_filtered = dedup_umis_4.combine(raw_bams, by:4)


// process getReadDistributionDifference{
//
//   storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/read_distribution_nvsf"
//
//   input:
//   set complete_id, dataset_name, dayPostInfection, tissue, sample,
//       file(umi_bam),
//       file(umi_bai),
//       dataset_name_2, dayPostInfection_2, tissue_2, sample_2,
//           file(bam),
//           file(bai) from raw_and_filtered
//   file(gene_annotation) from gene_annotation_rheMac_bed_2.collect()
//
//   output:
//   file "*" into read_distribution_channel_diff
//
//   script:
//   """
//   samtools view -h -o umi.sam ${umi_bam}
//   samtools view -h -o raw.sam  ${bam}
//   samtools sort  umi.sam > sorted_umi.sam
//   samtools sort  raw.sam > sorted_raw.sam
//   bam diff --in1 sorted_raw.sam --in2 sorted_umi.sam --out diff.sam
//   read_distribution.py -i diff_only1_sorted_raw.sam -r ${gene_annotation} > ${complete_id}.UMI.f3.q60.DIFF.read_distribution.txt
//   """
// }
// process geneAbundanceUmis{
//
//   storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/rnaseqc"
//
//   input:
//   set dataset_name, dayPostInfection, tissue, sample, complete_id,
//       file(bam),
//       file(bai) from dedup_umis_3
//   file annotation_gtf from gtf_transcriptonly_channel.collect()
//   file genome_fasta from reference_genome_ch.collect()
//   file genome_fai from reference_genome_fai_ch.collect()
//
//   output:
//   file("output/*") into output_rnaseqc_umis
//
//   script:
//   file_prefix = get_file_name_no_extension(bam.name)
//   """
//   mkdir output
//   # QC metrics
//   /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -Xmx4g -jar /opt/RNA-SeQC_1.1.9/RNA-SeQC.jar -n 1000 \
//               -s "${file_prefix}|${bam}|${complete_id}" -t ${annotation_gtf} \
//               -r ${genome_fasta} -o output/ -noDoC -strictMode -gatkFlags \
//               --allow_potentially_misencoded_quality_scores
//
//   """
// }
//
// process geneAbundanceMds{
//
//   storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/rnaseqc"
//
//   input:
//   set dataset_name, dayPostInfection, tissue, sample, complete_id,
//       file(bam),
//       file(bai) from md_bams_3
//   file annotation_gtf from gtf_transcriptonly_channel_2.collect()
//   file genome_fasta from reference_genome_ch_2.collect()
//   file genome_fai from reference_genome_fai_ch_2.collect()
//
//   output:
//   file("output/*") into output_rnaseqc_mds
//
//   script:
//   file_prefix = get_file_name_no_extension(bam.name)
//   """
//   mkdir output
//   # QC metrics
//   /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -Xmx4g -jar /opt/RNA-SeQC_1.1.9/RNA-SeQC.jar -n 1000 \
//               -s "${file_prefix}|${bam}|${complete_id}" -t ${annotation_gtf} \
//               -r ${genome_fasta} -o output/ -noDoC -strictMode -gatkFlags \
//               --allow_potentially_misencoded_quality_scores
//
//   """
// }
//
// process geneAbundanceFilter{
//
//   storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/rnaseqc"
//
//   input:
//   set dataset_name, dayPostInfection, tissue, sample, complete_id,
//       file(bam),
//       file(bai) from filtered_merged_bams_5
//   file annotation_gtf from gtf_transcriptonly_channel_3.collect()
//   file genome_fasta from reference_genome_ch_3.collect()
//   file genome_fai from reference_genome_fai_ch_3.collect()
//
//   output:
//   file("output/*") into output_rnaseqc_filter
//
//   script:
//   file_prefix = get_file_name_no_extension(bam.name)
//   """
//   mkdir output
//   # QC metrics
//   /usr/lib/jvm/java-1.7.0-openjdk-amd64/bin/java -Xmx4g -jar /opt/RNA-SeQC_1.1.9/RNA-SeQC.jar -n 1000 \
//               -s "${file_prefix}|${bam}|${complete_id}" -t ${annotation_gtf} \
//               -r ${genome_fasta} -o output/ -noDoC -strictMode -gatkFlags \
//               --allow_potentially_misencoded_quality_scores
//   """
// }


// =================================================================================
/*  # Gene-level expression, based on a collapsed annotation GTF
*  #java -Xmx6g -jar /apps/RNA-SeQC/1.1.9/RNA-SeQC.jar -n 1000 -s "${complete_id}|${bam}|${complete_id}" -t ${annotation_gtf} -r ${genome_fasta} -o output/ -noDoC -strictMode
*  DE NOVO TRANSCRIPTOME ASSEMBLY : StringTie/bin/bash -c "/src/run_rnaseqc.py ${bam} ${annotation_gtf} ${genome_fasta} ${complete_id} --output_dir output"
*/

params.transcriptome_assembly = true
process DeNovoAssembly{

  tag "${complete_id}"
  storeDir "${params.output_dir}/04_stringtie/$dataset_name/$tissue/$dayPostInfection/$sample"

  when:
  params.transcriptome_assembly

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(bam), file(bai) from dedup_umi_bams
  file gtf from gtf_rheMac_channel.collect()

  output:
  val dataset_name into dataset_name_ch
  file "*.stringtie.gtf" into stringTie_channel

  script:
  file_prefix = get_file_name_no_extension(bam.name)
  """
  stringtie ${bam} -G ${gtf} -o ${file_prefix}.stringtie.gtf -A -B --fr -f 0.1 -m 200 -a 10 -j 1 -c 0.1 -g 10 -p 2
  """
}

/*
* Merge all the assemblies Reference Guided
*/
process StringTie_Merge_Reference_Guided{

  cpus 1
  storeDir "${params.output_dir}/04_stringtie/$dataset_name"

  when:
  params.transcriptome_assembly

  input:
  val dataset_name from dataset_name_ch
  file(stringtie_gtfs) from stringTie_channel.collect()
  file reference_gtf from gtf_rheMac_channel2

  output:
  file "${dataset_name}_stringtie_merged_reference_guided.gtf" into (merged_denovo_assmebly, merged_de_novo_assembly_2)
  val dataset_name into dataset_name_ch2

  script:
  """
  stringtie --merge -p ${task.cpus} -o ${dataset_name}_stringtie_merged_reference_guided.gtf -G ${reference_gtf} ${stringtie_gtfs}
  """
}


process gffCompare{

  storeDir "${params.output_dir}/04_stringtie/$dataset_name/01_gffCompare"

  input:
  file merged_gtf from merged_denovo_assmebly
  file reference_gtf from gtf_rheMac_channel3
  val dataset_name from dataset_name_ch2

  output:
  file("merged*") into gff_compare_output_channel

  script:
  """
  gffcompare -r ${reference_gtf} -G -o merged ${merged_gtf}
  """
}


merged_de_novo_assembly_2.into{ merged_de_novo_assembly_3; merged_de_novo_assembly_4; merged_de_novo_assembly_5}


process StringTie_abundances_umis{

  cpus 8
  storeDir "${params.output_dir}/04_stringtie_counts/$dataset_name/$tissue/$dayPostInfection/$sample/umis"

  input:
  file stringtie_merged from merged_de_novo_assembly_3.collect()
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bam),
      file(bai) from dedup_umis_3

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,file_prefix, file("${file_prefix}_abundance.gtf") into abundances_stringtie_umis
  file "*.ctab" into ctabs_umis
  script:
  file_prefix = get_file_name_no_extension(bam.name)
  """
  stringtie -e -B -p ${task.cpus} -G ${stringtie_merged} -o ${file_prefix}_abundance.gtf ${bam}
  """
}

process StringTie_abundances_mds{

  cpus 8
  storeDir "${params.output_dir}/04_stringtie_counts/$dataset_name/$tissue/$dayPostInfection/$sample/mds"

  input:
  file stringtie_merged from merged_de_novo_assembly_4.collect()
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bam),
      file(bai) from md_bams_3

  output:
  file "${file_prefix}_abundance.gtf" into abundances_stringtie_mds
  file "*.ctab" into ctabs_mds

  script:
  file_prefix = get_file_name_no_extension(bam.name)
  """
  stringtie -e -B -p ${task.cpus} -G ${stringtie_merged} -o ${file_prefix}_abundance.gtf ${bam}
  """
}

process StringTie_abundances_filter{

  cpus 8
  storeDir "${params.output_dir}/04_stringtie_counts/$dataset_name/$tissue/$dayPostInfection/$sample/filter"

  input:
  file stringtie_merged from merged_de_novo_assembly_5.collect()
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bam),
      file(bai) from filtered_merged_bams_5

  output:
  file "${file_prefix}_abundance.gtf" into abundances_stringtie_filter
  file "*.ctab" into ctabs_filter

  script:
  file_prefix = get_file_name_no_extension(bam.name)
  """
  stringtie -e -B -p ${task.cpus} -G ${stringtie_merged} -o ${file_prefix}_abundance.gtf ${bam}
  """
}


process htseqCount_after_strintie{

  tag "${complete_id}"
  storeDir "${params.output_dir}/06_quantification/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,file(md_bam), file(md_metrics_bam) from marked_d_bams
  file gtf from gtfChannel.collect()

  output:
  file "${md_file_prefix}.HTseq.gene_counts.tab" into htseqCountsMd_channel

  script:
  md_file_prefix = get_file_name_no_extension(md_bam.name)
  """
  # Calc the counts for the md duplicates
  htseq-count -f bam -r pos -s yes -t gene -i gene_id ${md_bam} ${gtf} > ${md_file_prefix}.HTseq.gene_counts.tab
  """
}

params.script_prepde="${baseDir}/scripts/prepDE_mod.py"
Channel.fromPath("${params.script_prepde}").into{ prepdescript }

process prepareCountMatrices{

  storeDir "${params.output_dir}/04_StringtieAbundances/$dataset_name/$tissue/$dayPostInfection/$sample/counts"


  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file_prefix, file(abundance_gtf) from abundances_stringtie_umis
  file script from prepdescript.collect()

  output:
  file "*" into out_prepde_umis

  script:
  """
  python2 ${script} -i ${abundance_gtf} --sampleid ${complete_id} -g "${file_prefix}_gene_count_matrix.csv" -t "${file_prefix}_transcript_count_matrix.csv"
  """
}

workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

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
