params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/"


params.known_mrna = "${params.output_dir}/01_PreliminaryFiles_rheMac10/Macaca_mulatta.Mmul_10.100_coding_rna.fa"
params.known_lncrna = "${params.output_dir}/01_PreliminaryFiles_rheMac10/Macaca_mulatta.Mmul_10.100_noncoding_rna.fa"
Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }

reference_genome = Channel.fromPath("/gpfs/projects/bsc83/Data/assemblies/ensembl/release-100/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa").collect()
reference_genome_index = Channel.fromPath("/gpfs/projects/bsc83/Data/assemblies/ensembl/release-100/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa.fai").collect()

reference_gtf = Channel.fromPath("/gpfs/projects/bsc83/Data/gene_annotation/ensembl_release100/rheMac10/Macaca_mulatta.Mmul_10.100.gtf").collect()

inputgtf = Channel.fromPath("${params.output_dir}/01_PreliminaryFiles_rheMac10/Macaca_mulatta.Mmul_10.100_known_lncrna.gtf").collect()

log.info "=============================================="
log.info "            lncRNA annotation"
log.info "=============================================="

process feelnc_filter{

  storeDir "${params.output_dir}/00_feelNC_benchmark_macaque_filter/Input_known_lncrnas/"

  input:
  file gtf from inputgtf
  file reference_gtf from reference_gtf

  output:
  file("candidate_lncRNA.gtf") into candidates

  script:
  """
  FEELnc_filter.pl -i ${gtf} \
                   -a ${reference_gtf} \
                   -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
  """
}

process feelnc_codpot{

  cpus 16
  storeDir "${params.output_dir}/00_feelNC_benchmark_macaque_filter/Test_pl/"

  input:
  file candidate_lncrna from candidates
  file known_mrna from known_mrna_channel
  file known_lncrna from known_lncrna_channel
  file reference_genome
  file reference_genome_index

  output:
  file("*") into coding_potentials

  script:
  """
  FEELnc_codpot.pl -i ${candidate_lncrna} -a ${known_mrna} -l ${known_lncrna} \
                   -g ${reference_genome}  --proc ${task.cpus}
  """

}
