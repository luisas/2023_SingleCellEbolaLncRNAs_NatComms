
//params.data_folder = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data"
params.proj_data_folder = "/gpfs/projects/bsc83/Data"
params.dir_ebola_data = "${params.proj_data_folder}/Ebola"
params.data_folder = "${params.dir_ebola_data}/01_Ebola-RNASeq"
params.assembly_name = "rheMac10_EBOV-Kikwit"
params.output_dir = "${params.dir_ebola_data}/01_Ebola-RNASeq/02_RNA-Seq_rheMac10/"


params.ref_gtf = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/${params.assembly_name}.gtf"
params.merged_gtf = "${params.data_folder}/02_RNA-Seq_rheMac10/04_stringtie/stringtie_merged_reference_guided.gtf"

// TRAINING DATA
// Trained with human data!
params.known_mrna = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownMrnas_nochr.gtf"
params.known_lncrna = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownlncRNAs_nochr.gtf"

params.rhesus_genome = "${params.dir_ebola_data}/01_Ebola-RNASeq/01_PreliminaryFiles_rheMac10/reference_assembly/rheMac10_EBOV-Kikwit.fa"

Channel.fromPath("${params.rhesus_genome}").into{ fasta_reference_channel; fasta_reference_channel2}
Channel.fromPath("${params.merged_gtf}").into{ merged_gtf_channel; merged_gtf_channel_2 }
Channel.fromPath("${params.ref_gtf}").into{ ref_gtf_channel; ref_gtf_channel_2 }
Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }


log.info "=============================================="
log.info "            lncRNA annotation"
log.info "=============================================="

process feelnc_filter{

  storeDir "${params.output_dir}/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc"

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


process gffCompare{

  storeDir "${params.output_dir}/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc/01_gffcompare"

  input:
  set index, codpot from coding_potentials
  file reference_gtf from ref_gtf_channel_2

  output:
  file("merged*") into gff_compare_output_channel

  script:
  """
  gffcompare -R -r ${reference_gtf} -o merged ${codpot}/candidate_lncRNA.gtf.lncRNA.gtf
  """
}

// process feelnc_classifier{
//   storeDir "${params.output_dir}/05_lncrnaAnnotation/feelnc_gencode_linc"
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
//   perl /apps/FEELNC/0.1.1/scripts/FEELnc_classifier.pl -i feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf -a  ${reference_gtf} -l log.txt > lncRNA_classes.txt
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
