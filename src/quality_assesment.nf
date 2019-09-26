#!/usr/bin/env nextflow

/*
 * Nextflow pipeline for quality assessment of RNA-seq data.
 */
 log.info "=============================================="
 log.info "RNASeq analysis pipeline for Zyagen Samples"
 log.info "=============================================="

// ------------  INPUT PARAMETERS -------------
//params.dataset_bam_dir = "/gpfs/projects/bsc83/Ebola/00_RawData/pardis_shared_data/sabeti-txnomics/alin/190713_Zyagen-longRNA/tmp/00_demux/bams_per_lane/*/"
params.dataset_bam_dir = "/gpfs/scratch/bsc83/bsc83024/test_dataset/bams_per_lane/*/"

// Folder where the output directories of the pipeline will be placed
//params.output_dir = "/gpfs/projects/bsc83/Ebola/data/02_RNA-Seq/"
params.output_dir = "/gpfs/projects/bsc83/Ebola/test_data/02_RNA-Seq_2/"

// If set to true the quality assessment will be computed with fastqc
params.fastqc=true
params.mapping_statistics=true
// params needed for mapping
params.reference_assembly = "/gpfs/projects/bsc83/Ebola/data/01_PreliminaryFiles/reference_assembly/rheMac8_EBOV-Kikwit.fa"
// TODO change naming
params.preliminary_files_dir="/gpfs/projects/bsc83/Ebola/data/01_PreliminaryFiles"
params.assembly = "${params.preliminary_files_dir}/indexes/hisat2/rheMac8_EBOV-Kikwit"
params.assembly_prefix = "${params.assembly}".tokenize('/')[-1]
// splicing
params.ss = "${params.preliminary_files_dir}/gene_annotations/rheMac8_EBOV-Kikwit.ss.txt"
// Gene annotation
params.gene_annotation = "${params.preliminary_files_dir}/gene_annotations/rheMac8_EBOV-Kikwit.bed"

params.dict = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_prefix}.dict"
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
                           it.baseName.split('_')[0] + "_" + it.baseName.split('_')[1] +"_"+ it.baseName.split('_')[2].split('-')[0] +"_"+ it.baseName.split('_')[2].split('-')[1] +"_"+ "l" +  it.parent.name.split('\\.')[1],
                           it ) }

unmapped_bams = Channel
                .fromPath("${params.dataset_bam_dir}/*_long.bam")
                .ifEmpty('bam files directory is empty')
                .map{ tuple(it.baseName.split('_')[0] + "_" + it.baseName.split('_')[1] +"_"+ it.baseName.split('_')[2].split('-')[0] +"_"+ it.baseName.split('_')[2].split('-')[1] +"_"+ "l" +  it.parent.name.split('\\.')[1],
                            it)}

// Channel for the Assembly
assemblyForMapping = Channel.fromPath("${params.assembly}*")
reference_assembly_channel = Channel.fromPath("${params.reference_assembly}*")
// Channel for Annotation
annotationForMapping = Channel.fromPath("${params.ss}")

gene_annotation_channel = Channel.fromPath("${params.gene_annotation}")
dictionary_channel = Channel.fromPath("${params.dict}")


/*
* STEP 1: Convert the unmapped bam files into fastq files.
* Raw reads are stored in the bam format but are not mapped yet. They need
* to be converted into fastq format to proceed with the pipeline.
*/
process convert_bam_to_fastq {

    tag "${complete_id}"
    //publishDir  "${params.output_dir}/01_fastq/$dataset_name/$tissue/$dayPostInfection/$sample"
    storeDir "${params.output_dir}/01_fastq/${dataset_name}/$tissue/$dayPostInfection/$sample"

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


/*
* STEP 2.2: Map, sort and index
*/
process mapping_hisat{

  label 'big_mem'
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  file assembly from assemblyForMapping.collect()
  file ss from annotationForMapping.collect()
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
  hisat2-align-s -x ${params.assembly_prefix} -1  ${fastq_1} -2 ${fastq_2} --known-splicesite-infile ${ss} --novel-splicesite-outfile ${complete_id}.novel_ss.txt --downstream-transcriptome-assembly  --time --summary-file ${complete_id}.hisat2_summary.txt --rna-strandness FR > ${complete_id}.sam
  """

}

/*
* STEP 2.2: sort and index
*/
process sort_and_index_bam{

  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"


  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(ss), file(summary),
      file(sam) from mapped_sam

  output:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(ss), file(summary),
      file(sam), file("${complete_id}.bam") into (sorted_and_index_bam, sorted_indexed_bams_for_stats)


  script:
  """
  java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar SortSam INPUT=${sam} OUTPUT=${complete_id}.bam SORT_ORDER=coordinate
  """


}

/*
* Add UMIs
*/
//samtools sort -T $TMPDIR/${complete_id}.tmp -n -O bam -o ${complete_id}.bam ${sam}
//samtools index ${complete_id}.bam
// getting by complete_id
// See how channel were built - complete id is always in the -2 position
unmapped_and_mapped_bams = sorted_and_index_bam.join(unmapped_bams, by:-2)

process MergeBamAlignment{

  tag "${complete_id}"
  storeDir "${params.output_dir}/04_hisatTest/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(ss), file(summary),
      file(sam), file(bam),
      complete_id_unmapped, file(unmapped_bam) from unmapped_and_mapped_bams
  file assembly from reference_assembly_channel
  file dict from dictionary_channel

  output:
  file ("${complete_id}.UMI.bam") into merged_bam_alignments_channel

  script:
  """
  java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar MergeBamAlignment UNMAPPED=${unmapped_bam} ALIGNED=${bam} O=${complete_id}.UMI.bam R=${assembly} SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 ORIENTATIONS=FR CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT
  """
}


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

// Merge lanes in one
process mergeBAMperLane{

  tag "${complete_id}"
  storeDir "${params.output_dir}/04_hisatTest/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(samples) from groups_lanes

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file("${complete_id}.UMI.bam"), file("${complete_id}.UMI.bai") into merged_bylanes

  script:
  """
  #samtools merge -f ${key}.UMI.bam  ${outpath}/${sample_name}.l1.UMI.bam  ${outpath}/${sample_name}.l2.UMI.bam  ${outpath}/${sample_name}.l3.UMI.bam  ${outpath}/${sample_name}.l4.UMI.bam
  samtools merge -f ${complete_id}.UMI.bam  ${samples}
  # Merge mapped bam files with UMI info
  samtools index  ${complete_id}.UMI.bam
  """
}


// Filter put bams - now it is stringent

process filterBam{

  tag "${complete_id}"
  storeDir "${params.output_dir}/04_hisatTest/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(merged_bam), file(bai_mereged_bam) from merged_bylanes

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file("${complete_id}.UMI.f3.q60.bam") into (filtered_merged_bams_1, filtered_merged_bams_2)

  script:
  """
  samtools view -b -f3 -q 60  ${merged_bam} > ${complete_id}.UMI.f3.q60.bam
  samtools index ${complete_id}.UMI.f3.q60.bam
  """
}


// Standard marking duplicates -- prob we need to delete it

process MarkDuplicates{

  tag "${complete_id}"
  storeDir "${params.output_dir}/04_hisatTest/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(filtered_bam) from filtered_merged_bams_1

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,file("${complete_id}.md.bam"), file("${complete_id}.md_metrics.bam") into marked_d_bams

  script:
  """
  java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar MarkDuplicates I=${filtered_bam} O=${complete_id}.md.bam M=${complete_id}.md_metrics.bam REMOVE_DUPLICATES=false ASSUME_SORTED=true  CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
  """
}


//dedup duplicates using UMI

process dedupUmi{

  tag "${complete_id}"
  storeDir "${params.output_dir}/04_hisatTest/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(filtered_bam) from filtered_merged_bams_2


  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file("${complete_id}.UMI.f3.q60.umi_dedup.bam"),
      file("${complete_id}.UMI.f3.q60.umi_dedup.bai") into dedup_umi_bams

  script:
  """
  umi_tools dedup --stdin=${filtered_bam} --log=${complete_id}.umi_tools.dedup.log --output-stats=${complete_id} --extract-umi-method=tag --umi-tag=RX --paired --method directional  > ${complete_id}.UMI.f3.q60.umi_dedup.bam
  samtools index ${complete_id}.UMI.f3.q60.umi_dedup.bam
  """

}


process getCountsUMIs{

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(umi_bam),
      file(umi_bai) from dedup_umi_bams

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

  echo -e "Total\t"\${n} > counts/${complete_id}.n_reads.txt
  echo -e "Autosomes\t"\${n_autosomes} >> counts/${complete_id}.n_reads.txt
  echo -e "chrX\t"\${n_chrX} >> counts/${complete_id}.n_reads.txt
  echo -e "chrY\t"\${n_chrY} >> counts/${complete_id}.n_reads.txt
  echo -e "chrM\t"\${n_chrM} >> counts/${complete_id}.n_reads.txt
  echo -e "EBOV\t"\${n_EBOV} >> counts/${complete_id}.n_reads.txt
  """
}

/*
* STEP 3: stats
*/

// process mapping_stats{
//
//   tag "${complete_id}"
//   storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"
//
//   when:
//   params.mapping_statistics
//
//   input:
//   file gene_annotation from gene_annotation_channel.collect()
//   set lane,lane_number, dataset_name, dayPostInfection, tissue, sample, complete_id,
//       file(ss), file(summary),
//       file(sam), file(bam) from sorted_indexed_bams_for_stats
//
//   output:
//   set file("infer_experiment/*"), file("clipping_profile/*") into reseqstats
//
//
//   script:
//   """
//
//   geneBody_coverage.py -i ${bam_files} -r /gpfs/projects/bsc83/Ebola/data/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8_EBOV-Kikwit/RefSeq/rheMac8_EBOV-Kikwit.bed -o geneBody_coverage
//
//   ## INFER EXPERIMENT
//   mkdir infer_experiment
//   infer_experiment.py -r ${gene_annotation}  -s 1000000 -i ${bam} > infer_experiment/${complete_id}.inferred_strandness.txt
//
//
//   # CLIPPING PROFILE
//   # hisat2 only 3 mapping qualities: 0 (unmapped), 1(multiply mapped) and 60(uniquely mapped)
//   mkdir clipping_profile
//   clipping_profile.py -i ${bam} -o clipping_profile/${complete_id} -q 60 -s PE
//
//   #FPKM count
//   mkdir fpkmCount
//   FPKM_count.py -i ${bam}-o fpkmCount/${complete_id} -r ${gene_annotation} -d '1++,1--,2+-,2-+' -u -e -q 60 --single-read=0
//
//
//   # MARK duplicates
//   # duplicates are marked with flagstat 1024
//   mdkir markedDuplicates
//   java -Xmx10g -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar MarkDuplicates I=${bam} O=markedDuplicates/${complete_id}.md.bam M=markedDuplicates/${complete_id}.md_metrics.bam REMOVE_DUPLICATES=false ASSUME_SORTED=true  CREATE_INDEX=true MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
//   """
//
// }

def remove_lane_from_id(String id){
  return id.split("_").init().join("_")
}
