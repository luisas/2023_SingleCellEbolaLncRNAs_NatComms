



log.info "=============================================="
log.info " Merge stringtie assemblies and perform prediction "
log.info "=============================================="


// BaseFolders
params.prefix = "rheMac10_EBOV-Kikwit_UCSC"
params.prefix_data = "/gpfs/projects/bsc83/Data"
params.output_dir_preliminary = "${params.prefix_data}/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/"
params.output_dir_name = "01_RNA-Seq_ribodepl"
params.output_dir = "${params.prefix_data}/Ebola/99_BroadAnnotation/"


// Reference annotation
params.reference_annotated = "${params.output_dir_preliminary}/gene_annotations/UCSC/rheMac10_EBOV-Kikwit_UCSC.gtf"
gtfChannel = Channel.fromPath("${params.reference_annotated}")
params.bed = "${params.output_dir_preliminary}/gene_annotations/UCSC/rheMac10_EBOV-Kikwit_UCSC.bed12"
gtfChannel_ucsc = Channel.fromPath("${params.output_dir_preliminary}/gene_annotations/rheMac10_EBOV-Kikwit_UCSC.gtf")

params.reference_fasta = "${params.output_dir_preliminary}/reference_assembly/rheMac10_EBOV-Kikwit_UCSC.fa"
ref_fasta = Channel.fromPath("${params.reference_fasta}")

gtfChannel.into{gtfChannel1; gtfChannel2; gtfChannel3; gtfChannel4;gtfChannel5;  }

// Stringtie Files
stringTie_channel = Channel.fromPath("${params.prefix_data}/Ebola/00_RawData/BroadTranscriptomesComplete/stringtie2-gtfs/stringtie2-gtfs/*.gtf")


params.known_mrna = "${params.output_dir_preliminary}/Homo_sapiens.GRCh38.83_antisense.lincRNA_learning5k.fa"
params.known_lncrna = "${params.output_dir_preliminary}/Homo_sapiens.GRCh38.83_protein_coding_learning5k.fa"

files = Channel.fromPath("/gpfs/projects/bsc83/Data/Ebola/00_RawData/BroadTranscriptomesComplete/bams/*")
                .ifEmpty("No bams found")
                .map { tuple(it.baseName, it) }


files.into{ bams; }

bed_channel = Channel.fromPath("${params.bed}")

Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }


scriptConversion =  Channel.fromPath("${baseDir}/scripts/convert_annotation.R").collect()


process runRSeQC{

  storeDir "${params.output_dir}/03b_rseqc/$dataset_name/$tissue/$dayPostInfection/$sample/rseqc"

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
