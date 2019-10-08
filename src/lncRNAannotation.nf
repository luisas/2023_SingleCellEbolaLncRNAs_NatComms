
params.data_folder = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data"
params.assembly_name = "rheMac8_EBOV-Kikwit"

params.output_dir = "/gpfs/projects/bsc83/Projects/Ebola/data/02_RNA-Seq/"
params.ref_gtf = "${params.data_folder}/02_RNA-Seq/01_PreliminaryFiles/gene_annotations/${params.assembly_name}.gtf"
params.merged_gtf = "${params.data_folder}/02_RNA-Seq/04_stringtie/stringtie_merged_reference_guided.gtf"

params.known_mrna = "${params.data_folder}/02_RNA-Seq/01_PreliminaryFiles/rheMac8_EBOV-Kikwit_nolong.gtf"
params.known_lncrna = "${params.data_folder}/02_RNA-Seq/01_PreliminaryFiles/rheMac8.ensembl_release97_knownlncrna.gtf"

params.ref_fasta_mac = "/gpfs/projects/bsc83/Data/assemblies/ensembl/release-96/Mmul_8.0.1/Macaca_mulatta.Mmul_8.0.1.dna.toplevel.fa"

Channel.fromPath("${params.ref_fasta_mac}").into{ fasta_reference_channel}
Channel.fromPath("${params.merged_gtf}").into{ merged_gtf_channel; merged_gtf_channel_2 }
Channel.fromPath("${params.ref_gtf}").into{ ref_gtf_channel }
Channel.fromPath("${params.known_mrna}").into{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").into{ known_lncrna_channel }


log.info "=============================================="
log.info "            lncRNA annotation"
log.info "=============================================="

process feelnc_filter{

  storeDir "${params.output_dir}/05_lncrnaAnnotation/feelnc"

  input:
  file merged_gtf from merged_gtf_channel
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

  storeDir "${params.output_dir}/05_lncrnaAnnotation/feelnc"

  input:
  file candidate_lncrna from candidates
  file known_mrna from known_mrna_channel
  file known_lncrna from known_lncrna_channel

  output:
  file("*") into coding_potentials

  script:
  """
  FEELnc_codpot.pl -i ${candidate_lncrna} -a ${known_mrna} -l ${known_lncrna}
  """

}

// --------------
// lncADeep

process get_fasta_from_gtf{

  storeDir "${params.output_dir}/05_lncrnaAnnotation/lncadeep"

  input:
  file gtf from merged_gtf_channel_2
  file fasta_reference from fasta_reference_channel

  output:
  file("extracted.fa") into fasta_to_predict_channel

  script:
  """
  bedtools getfasta -fi ${fasta_reference} -bed ${gtf} -fo extracted.fa
  """
}

process lncadeep{

  cpu 48
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
