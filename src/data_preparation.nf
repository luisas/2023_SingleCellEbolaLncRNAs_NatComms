// Nextflow pipeline for data preparation  .
// possibly to be later merged with the analysis pipeline

params.prefix_rawdata = "/gpfs/projects/bsc83/Ebola/00_RawData/"
params.rhesus_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8/rheMac8.fa"
params.ebov_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.fa"

rhesus_genome_channel = Channel.fromPath("${params.rhesus_genome}")
ebov_genome_channel = Channel.fromPath("${params.ebov_genome}")


params.output_dir = "/gpfs/scratch/bsc83/bsc83024/test_pipeline"
/*
* Merge Assemblies (+ change name to match with annotation in ebov_genome)
* GTF gene annotation EBOV annotates as EBOV_Kikwit while in fasta file as
* Replace EBOV_Kikwit for chrom name in fasta file so it matches with annotation
*/
process merge_assemblies {

    storeDir "${params.output_dir}/reference_assembly"

    input:
    file rheMac8 from rhesus_genome_channel
    file ebov from ebov_genome_channel

    output:
    file "rheMac8_EBOV-Kikwit.fa" into merged_assembly

    script:
    """
    cat ${rheMac8} > rheMac8_EBOV-Kikwit.fa
    sed 's/KU182905.1/EBOV_Kikwit/g' ${ebov} >> rheMac8_EBOV-Kikwit.fa
    """

}
