



log.info "=============================================="
log.info " Merge stringtie assemblies and perform prediction "
log.info "=============================================="

// Directories
params.prefix_data = "/gpfs/projects/bsc83/Data"
params.output_dir = "${params.prefix_data}/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_ribodepl"


//Bed file
params.bed = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/gene_annotations/rheMac10_EBOV_and_novel_genenames.bed12"
bed_channel = Channel.fromPath("${params.bed}").collect()


// Bam files
bams = Channel.fromPath("/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/02_RNA-Seq_ribodepl/03_hisat/*/*/*/*/*.umi_dedup.bam")
                .ifEmpty("No bams found")
                .map { tuple(it.baseName, it) }



process runRSeQC{

  cpus 8
  storeDir "${params.output_dir}/03c_rseqc/${complete_id}"

  input:
  set complete_id,
      file(bam) from bams
  file(bed12) from  bed_channel.collect()

  output:
  file "*" into read_distribution_channel

  script:
  """
  read_distribution.py -i ${bam} -r ${bed12} > ${complete_id}.UMI.f3.q60.read_distribution.txt
  infer_experiment.py -i ${bam} -r ${bed12} > ${complete_id}.infer_experiment.txt
  junction_annotation.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12}
  bam_stat.py -i ${bam} > ${complete_id}.bam_stat.txt
  junction_saturation.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12} > ${complete_id}.junction_annotation_log.txt
  inner_distance.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12}
  read_duplication.py -i ${bam} -o ${complete_id}.read_duplication
  """
}
