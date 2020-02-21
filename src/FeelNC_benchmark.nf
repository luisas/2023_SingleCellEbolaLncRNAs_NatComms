
// FeelNC is run on top of known 

params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.output_dir_name = "02_RNA-Seq_external"
params.output_dir = "${params.dirData}/01_Ebola-RNASeq_all/${params.output_dir_name}/"

params.data_folder = "${params.dirData}/01_Ebola-RNASeq"
params.known_mrna = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownMrnas_nochr.gtf"
params.known_lncrna = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownlncRNAs_nochr.gtf"
Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }


Channel.fromPath("${params.dirData}").set{fasta_reference_channel}

log.info "=============================================="
log.info "            lncRNA annotation"
log.info "=============================================="

process feelnc_filter{

  storeDir "${params.output_dir}/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc"

  input:
  file merged_gtf from merged_de_novo_assembly_3
  file reference_gtf from ref_gtf_channel

  output:
  file("candidate_lncRNA.gtf") into candidates

  script:
  """
  FEELnc_filter.pl -i ${merged_gtf} \
                   -a ${reference_gtf} \
                   -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
  """
}

process feelnc_codpot{
  cpus 48
  storeDir "${params.output_dir}/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc"

  input:
  file candidate_lncrna from candidates
  file known_mrna from known_mrna_channel
  file known_lncrna from known_lncrna_channel
  file reference_genome from fasta_reference_channel

  output:
  set file("rheMac10_EBOV-Kikwit.fa.index"), file("feelnc_codpot_out") into coding_potentials

  script:
  """
  FEELnc_codpot.pl -i ${candidate_lncrna} -a ${known_mrna} -l ${known_lncrna} \
                   -g ${reference_genome} --proc ${task.cpus}
  """

}
