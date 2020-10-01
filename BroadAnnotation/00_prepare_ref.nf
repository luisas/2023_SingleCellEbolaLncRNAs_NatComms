
// BaseFolders
params.prefix = "rheMac10_EBOV-Kikwit_UCSC"
params.prefix_data = "/gpfs/projects/bsc83/Data"
params.output_dir_preliminary = "${params.prefix_data}/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/"


// Reference annotation
params.reference_annotated = "${params.output_dir_preliminary}/gene_annotations/rheMac10_EBOV-Kikwit.gtf"
gtfChannel = Channel.fromPath("${params.reference_annotated}")


params.reference_fasta = "${params.output_dir_preliminary}/reference_assembly/rheMac10_EBOV-Kikwit.fa"
ref_fasta = Channel.fromPath("${params.reference_fasta}")

/*
* Merge Assemblies of macaque and EBOV Virus to generate one merged assembly.
*/
process merge_assemblies {

    storeDir "${params.output_dir_preliminary}/reference_assembly"

    input:
    file rheMac from rhesus_genome_channel
    file ebov from ebov_genome_channel

    output:
    set file("${params.prefix}.fa"), file("${params.prefix}.fa.fai") into (merged_assembly,fasta_reference_channel,  reference_assembly_channel,  merged_assembly_for_dictionary )

    script:
    """
    cat ${rheMac} > ${params.prefix}.fa
    samtools faidx ${params.prefix}.fa
    """

}

/*
* Concatenate ensembl release annotation with EBOV gene annotation
*/
process merge_annotations{
  storeDir "${params.output_dir_preliminary}/gene_annotations"

  input:
  file ebov_annotation from ebov_annotation_channel
  file(rheMac_annotation) from rheMac_annotation_channel

  output:
  file("${params.prefix}.gtf") into (merged_annotation_channel,  merged_annotation_ch, merged_annotation_ch_2  )

  script:
  """
  cat  ${rheMac_annotation} ${ebov_annotation} > ${params.prefix}.gtf
  """

}
