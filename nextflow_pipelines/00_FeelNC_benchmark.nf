
// Filter is performed on the correct RheMac gtf
params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/"


params.known_mrna = "${params.output_dir}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownMrnas_nochr.gtf"
params.known_lncrna = "${params.output_dir}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownlncRNAs_nochr.gtf"
Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }

reference_genome = Channel.fromPath("/gpfs/projects/bsc83/Data/assemblies/ensembl/release-98/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa").collect()
reference_genome_index = Channel.fromPath("/gpfs/projects/bsc83/Data/assemblies/ensembl/release-98/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa.fai").collect()

reference_gtf = Channel.fromPath("/gpfs/projects/bsc83/Data/gene_annotation/ensembl_release98/rheMac10/Macaca_mulatta.Mmul_10.98.gtf").collect()

log.info "=============================================="
log.info "            lncRNA annotation"
log.info "=============================================="

process feelnc_filter{

  storeDir "${params.output_dir}/00_feelNC_benchmark"

  input:
  file reference_gtf

  output:
  file("candidate_lncRNA.gtf") into candidates

  script:
  """
  FEELnc_filter.pl -i ${reference_gtf} \
                   -a ${reference_gtf} \
                   -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
  """
}

process feelnc_codpot{
  cpus 48
  storeDir "${params.output_dir}/00_feelNC_benchmark"

  input:
  file candidate_lncrna from candidates
  file known_mrna from known_mrna_channel
  file known_lncrna from known_lncrna_channel
  file reference_genome
  file reference_genome_index

  output:
  set file("rheMac10_EBOV-Kikwit.fa.index"), file("feelnc_codpot_out") into coding_potentials

  script:
  """
  FEELnc_codpot.pl -i ${candidate_lncrna} -a ${known_mrna} -l ${known_lncrna} \
                   -g ${reference_genome} --proc ${task.cpus}
  """

}
