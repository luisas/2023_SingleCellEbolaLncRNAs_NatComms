
// BaseFolders
params.prefix = "rheMac10_EBOV-Kikwit_UCSC"
params.prefix_data = "/gpfs/projects/bsc83/Data"
params.output_dir_preliminary = "${params.prefix_data}/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/"

// Reference Annotation and Assembly - Macaque
params.rhesus_gtf = "${params.prefix_data}/gene_annotation/UCSC/rheMac10/rheMac10.gtf"
params.rhesus_genome = "${params.prefix_data}/assemblies/UCSC/rheMac10/rheMac10.fa"

// Ebola virus annotation and assembly
params.prefix_rawdata = "${params.prefix_data}/Ebola/00_RawData"
params.ebov_genome = "${params.prefix_rawdata}/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.fa"
params.ebov_gtf = "${params.prefix_rawdata}/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.gtf"


rhesus_genome_channel = Channel
                        .fromPath("${params.rhesus_genome}")
ebov_genome_channel = Channel
                      .fromPath("${params.ebov_genome}")
scripts=file("${params.scripts}")
rheMac_annotation_channel = Channel
                                  .fromPath("${params.rhesus_gtf}")
ebov_annotation_channel = Channel.fromPath("${params.ebov_gtf}")


params.scripts="${baseDir}/scripts/"

gtfToBed_script_ch = Channel
                      .fromPath("${params.scripts}/gtf2bed")
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
  storeDir "${params.output_dir_preliminary}/gene_annotations/UCSC/"

  input:
  file ebov_annotation from ebov_annotation_channel
  file(rheMac_annotation) from rheMac_annotation_channel

  output:
  file("${params.prefix}.gtf") into (merged_annotation_channel,  merged_annotation_ch, merged_annotation_ch_2  )

  script:
  """
  cat  ${rheMac_annotation} ${ebov_annotation} > ${params.prefix}_old.gtf
  sed 's/EBOV_Kikwit/KU182905.1/g' ${params.prefix}_old.gtf >> ${params.prefix}.gtf
  """

}


/*
* Convert GTF to BED
*/
process convert_gtf_to_bed12{
  storeDir "${params.output_dir_preliminary}/gene_annotations"

  input:
  file merged_annotation from merged_annotation_ch
  file gtf2bed from gtfToBed_script_ch

  output:
  file("${params.prefix}.bed") into bed_channel

  script:
  """
  perl ./gtf2bed $merged_annotation > ${merged_annotation.baseName}.bed
  """
}
