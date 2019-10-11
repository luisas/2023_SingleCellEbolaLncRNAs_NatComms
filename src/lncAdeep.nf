
//params.data_folder = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data"
params.data_folder = "/gpfs/projects/bsc83/Projects/Ebola/data"
params.assembly_name = "rheMac8_EBOV-Kikwit"
params.prefix_rawdata="/gpfs/projects/bsc83/Data/Ebola/00_RawData/"

params.output_dir = "/gpfs/projects/bsc83/Projects/Ebola/data/02_RNA-Seq/"
params.ref_gtf = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/${params.assembly_name}.gtf"
params.merged_gtf = "${params.data_folder}/02_RNA-Seq/04_stringtie/Zyagen/Zyagen_stringtie_merged_reference_guided.gtf"

//params.known_mrna = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/rheMac8_EBOV-Kikwit_nolong.gtf"
//params.known_lncrna = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/rheMac8.ensembl_release97_knownlncrna.gtf"
// TRAINING DATA
params.known_mrna = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/gencode.v26.GRCh38.annotation_knownMrnas.gtf "
params.known_lncrna = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/gencode.v26.GRCh38.annotation_knownlncRNAs.gtf"

params.rhesus_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8/rheMac8.fa"

Channel.fromPath("${params.rhesus_genome}").into{ fasta_reference_channel; fasta_reference_channel2}
Channel.fromPath("${params.merged_gtf}").into{ merged_gtf_channel; merged_gtf_channel_2 }
Channel.fromPath("${params.ref_gtf}").set{ ref_gtf_channel }
Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }


log.info "=============================================="
log.info "       lncRNA annotation w/ lncADeep "
log.info "=============================================="


process get_fasta_from_gtf{

  storeDir "${params.output_dir}/05_lncrnaAnnotation/lncadeep"

  input:
  file gtf from merged_gtf_channel_2
  file fasta_reference from fasta_reference_channel2

  output:
  file("extracted.fa") into fasta_to_predict_channel

  script:
  """
  bedtools getfasta -fi ${fasta_reference} -bed ${gtf} -fo extracted.fa
  """
}

process lncadeep{

  cpus 48
  storeDir "${params.output_dir}/05_lncrnaAnnotation/lncadeep"

  input:
  file(fasta_to_predict) from fasta_to_predict_channel

  output:
  file("lncadeep*") into lncadeep_output_channel

  script:
  """
  LncADeep.py -MODE lncRNA -f ${fasta_to_predict} -o lncadeep -th ${task.cpus}
  """
}
