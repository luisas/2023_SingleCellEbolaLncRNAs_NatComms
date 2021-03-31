




log.info "=============================================="
log.info " merge stringtie assemblies and perform prediction "
log.info "=============================================="


// BaseFolders
params.prefix = "rheMac10_EBOV-Kikwit"
params.prefix_data = "/gpfs/projects/bsc83/Data"
params.output_dir_preliminary = "${params.prefix_data}/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/"
params.output_dir_name = "02_RNA-Seq_ribodepl"
params.output_dir = "${params.prefix_data}/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/${params.output_dir_name}/"


// Reference annotation
params.reference_annotated = "${params.output_dir_preliminary}/gene_annotations/rheMac10_EBOV-Kikwit.gtf"
gtfChannel = Channel.fromPath("${params.reference_annotated}")


params.reference_fasta = "${params.output_dir_preliminary}/reference_assembly/rheMac10_EBOV-Kikwit.fa"
ref_fasta = Channel.fromPath("${params.reference_fasta}")

gtfChannel.into{gtfChannel1; gtfChannel2; gtfChannel3; gtfChannel4;gtfChannel5;  }

// Stringtie Files
stringTie_channel = Channel.fromPath("${params.output_dir}/04_stringtie/*/*/*/*/*stringtie.gtf")


// Training data for Feelnc
//params.known_mrna = "${params.output_dir_preliminary}/Homo_sapiens.GRCh38.100_coding_rna.fa"
//params.known_lncrna = "${params.output_dir_preliminary}/Homo_sapiens.GRCh38.100_noncoding_rna.fa"


params.known_mrna = "${params.output_dir_preliminary}/Homo_sapiens.GRCh38.83_antisense.lincRNA_learning5k.fa"
params.known_lncrna = "${params.output_dir_preliminary}/Homo_sapiens.GRCh38.83_protein_coding_learning5k.fa"

Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }


/*
* Merge all the assemblies Reference Guided
*/
process StringTie_Merge_Reference_Guided{

  cpus 1
  storeDir "${params.output_dir}/04_stringtie/"


  input:
  file(stringtie_gtfs) from stringTie_channel.collect()
  file reference_gtf from gtfChannel1

  output:
  file "stringtie_merged_reference_guided.gtf" into (merged_denovo_assmebly, merged_de_novo_assembly_2)

  script:
  """
  stringtie --merge -p ${task.cpus} -o stringtie_merged_reference_guided.gtf -G ${reference_gtf} ${stringtie_gtfs}
  """
}


/*
* Compare stringtie output with reference GTF
*/
process gffCompare2{

  storeDir "${params.output_dir}/04_stringtie_gffcompare/01_gffCompare"

  input:
  file merged_gtf from merged_denovo_assmebly
  file reference_gtf from gtfChannel2

  output:
  file("merged*") into gff_compare_output_channel2

  script:
  """
  gffcompare -R -r ${reference_gtf} -o merged ${merged_gtf}
  """
}


merged_de_novo_assembly_2.into{ merged_de_novo_assembly_3; merged_de_novo_assembly_4; merged_de_novo_assembly_5}



log.info "=============================================="
log.info "            lncRNA Annotation w/ FeelNC "
log.info "=============================================="

/*
* Deplete transcripts that are: short (<200 nucleotides), monoexonic, overlapping protein coding genes
*/

process feelnc_filter{

  storeDir "${params.output_dir}/05_feelNC_5k/"

  input:
  file merged_gtf from merged_de_novo_assembly_3
  file reference_gtf from gtfChannel3

  output:
  file("candidate_lncRNA.gtf") into candidates

  script:
  """
  FEELnc_filter.pl -i ${merged_gtf} \
                   -a ${reference_gtf} \
                   -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
  """
}

/*
* Predict novel lncRNAs
*/

process feelnc_codpot{
  //cpus 48
  cpus 16
  storeDir "${params.output_dir}/05_feelNC_5k"

  input:
  file candidate_lncrna from candidates
  file known_mrna from known_mrna_channel
  file known_lncrna from known_lncrna_channel
  file reference_genome from ref_fasta

  output:
  set file("rheMac10_EBOV-Kikwit.fa.index"), file("feelnc_codpot_out") into coding_potentials2

  script:
  """
  FEELnc_codpot.pl -i ${candidate_lncrna} -a ${known_mrna} -l ${known_lncrna} \
                   -g ${reference_genome} --proc ${task.cpus}
  """

}

coding_potentials2.into{coding_potentials; coding_potentials1; }

process feelnc_classifier{
  storeDir "${params.output_dir}/05_feelNC_5k/02_feelnc_classifier"

  input:
  set index, codpot from coding_potentials1
  file reference_gtf from gtfChannel5

  output:
  file("rheMac10_EBOV-Kikwit_lncRNA_classes.txt") into classification

  script:
  """
  FEELnc_classifier.pl -i ${codpot}/candidate_lncRNA.gtf.lncRNA.gtf -a  $reference_gtf > rheMac10_EBOV-Kikwit_lncRNA_classes.txt

  """

}

/*
* Compare prediction with reference GTF
*/
process gffCompare{

  storeDir "${params.output_dir}/05_feelNC_5k/01_gffcompare"

  input:
  set index, codpot from coding_potentials
  file reference_gtf from gtfChannel4

  output:
  file("merged*") into gff_compare_output_channel

  script:
  """
  gffcompare -R -r ${reference_gtf} -o merged ${codpot}/candidate_lncRNA.gtf.lncRNA.gtf
  """
}



workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

/*   -------------------------------
*           Groovy Functions
*    -------------------------------
*/

def remove_lane_from_id(String id){
  return id.split("_").init().join("_")
}

def get_file_name_no_extension(String filename){
  return filename.split("\\.").init().join('.')
}
