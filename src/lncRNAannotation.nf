
//params.data_folder = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data"
params.data_folder = "/gpfs/projects/bsc83/Projects/Ebola/data"
params.assembly_name = "rheMac8_EBOV-Kikwit"
params.prefix_rawdata="/gpfs/projects/bsc83/Data/Ebola/00_RawData/"

params.output_dir = "/gpfs/projects/bsc83/Projects/Ebola/data/02_RNA-Seq/"
params.ref_gtf = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/${params.assembly_name}.gtf"
//params.ref_gtf = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/${params.assembly_name}_mrnas.gtf"
//params.merged_gtf = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/${params.assembly_name}.gtf"
params.merged_gtf = "${params.data_folder}/02_RNA-Seq/04_stringtie/Zyagen/Zyagen_stringtie_merged_reference_guided.gtf"

//params.known_mrna = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/rheMac8_EBOV-Kikwit_nolong.gtf"
//params.known_lncrna = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/rheMac8.ensembl_release97_knownlncrna.gtf"
// TRAINING DATA
params.known_mrna = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/gencode.v26.GRCh38.annotation_knownMrnas.gtf"
params.known_lncrna = "${params.data_folder}/01_PreliminaryFiles/gene_annotations/gencode.v26.GRCh38.annotation_knownlncRNAs.gtf"

params.rhesus_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8/rheMac8.fa"

Channel.fromPath("${params.rhesus_genome}").into{ fasta_reference_channel; fasta_reference_channel2}
Channel.fromPath("${params.merged_gtf}").into{ merged_gtf_channel; merged_gtf_channel_2 }
Channel.fromPath("${params.ref_gtf}").into{ ref_gtf_channel; ref_gtf_channel_2 }
Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }


log.info "=============================================="
log.info "            lncRNA annotation"
log.info "=============================================="

process feelnc_filter{

  storeDir "${params.output_dir}/05_lncrnaAnnotation/feelnc_gencode_two_cut_offs"

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
  cpus 48
  storeDir "${params.output_dir}/05_lncrnaAnnotation/feelnc_gencode_two_cut_offs_pred"

  input:
  file candidate_lncrna from candidates
  file known_mrna from known_mrna_channel
  file known_lncrna from known_lncrna_channel
  file reference_genome from fasta_reference_channel

  output:
  file("*") into coding_potentials

  script:
  """
  FEELnc_codpot.pl -i ${candidate_lncrna} -a ${known_mrna} -l ${known_lncrna} \
                   -g ${reference_genome} --proc ${task.cpus} \
                   --spethres 0.96,0.96
  """

}




// process feelnc_classifier{
//   storeDir "${params.output_dir}/05_lncrnaAnnotation/feelnc_gencode_lncrna"
//
//   input:
//   file(cod_pot_dir) from coding_potentials
//   file reference_gtf from ref_gtf_channel_2
//
//   output:
//   file("lncRNA_classes.txt") into classification_ch
//
//   script:
//   """
//   FEELnc_classifier.pl -i ${cod_pot_dir}/candidate_lncRNA.gtf.lncRNA.gtf -a  ${reference_gtf} > lncRNA_classes.txt
//   """
//
// }




// // --------------
// // lncADeep
//
// process get_fasta_from_gtf{
//
//   storeDir "${params.output_dir}/05_lncrnaAnnotation/lncadeep"
//
//   input:
//   file gtf from merged_gtf_channel_2
//   file fasta_reference from fasta_reference_channel2
//
//   output:
//   file("extracted.fa") into fasta_to_predict_channel
//
//   script:
//   """
//   bedtools getfasta -fi ${fasta_reference} -bed ${gtf} -fo extracted.fa
//   """
// }

// process lncadeep{
//
//   cpu 48
//   storeDir "${params.output_dir}/05_lncrnaAnnotation/lncadeep"
//
//   input:
//   file(fasta_to_predict) from fasta_to_predict_channel
//
//   output:
//   file("lncadeep*") into lncadeep_output_channel
//
//   script:
//   """
//   LncADeep.py -MODE lncRNA -f ${fasta_to_predict} -o lncadeep -th ${task.cpus}
//   """
// }
