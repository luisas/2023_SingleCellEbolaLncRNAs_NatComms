#!/usr/bin/env nextflow

/*
 * TODO COPYRIGHT
 */

 log.info "=============================================="
 log.info "RNASeq analysis pipeline for Zyagen Samples"
 log.info "=============================================="

 /* --------------------------------------
  *  INPUT PARAMETERS defaults
  * --------------------------------------
  */
// Unmapped bam files folder - RAW DATA
params.dataset_bam_dir = "/gpfs/projects/bsc83/Data/Ebola/00_RawData/pardis_shared_data/sabeti-txnomics/alin/190713_Zyagen-longRNA/tmp/00_demux/bams_per_lane/*/"
//params.dataset_bam_dir = "/gpfs/scratch/bsc83/bsc83024/test_dataset/bams_per_lane/*/"

// Folder where the output directories of the pipeline will be placed
params.output_dir = "/gpfs/projects/bsc83/Projects/Ebola/data/02_RNA-Seq/"
//params.output_dir = "/gpfs/projects/bsc83/Projects/Ebola/test_data/04_test/"


params.preliminary_files_dir="/gpfs/projects/bsc83/Projects/Ebola/data/01_PreliminaryFiles"
params.assembly_name = "rheMac8_EBOV-Kikwit"

// This ones should be automatically retrieved but can be changed in need
params.reference_assembly = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.fa"
// Dictionary
params.dict = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.dict"
// Indexes
params.hisat2_indexes = "${params.preliminary_files_dir}/indexes/hisat2/${params.assembly_name}"
// Known splice sites
params.ss = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.ss.txt"
// Gene annotation
params.gene_annotation = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.bed"
params.gtf = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.gtf"
params.gtf_rheMac = "${params.preliminary_files_dir}/gene_annotations/rheMac8.ensembl_release97.gtf"
params.bed_rheMac = "/gpfs/projects/bsc83/Data/Ebola/00_RawData/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8/Ensembl/rheMac8.Ensembl.bed"
// If set to true the quality assessment will be computed with fastqc
params.fastqc=true

// -----------------------------------------------

/* Channel for all the unmapped bam files present in any of the dataset_bam_dir
* subdirectories.
* In order, it maps:
* - lane
* - lane_number
* - dataset_name
* - dayPostInfection
* - tissue
* - sample
* - complete_id (${dataset_name}_${dayPostInfection}_${tissue}_${sample}_l${lane_number})
* - the file itself.
*/
dataset_bam = Channel
              .fromPath("${params.dataset_bam_dir}/*_long.bam")
              .ifEmpty('bam files directory is empty')
              .map { tuple(it.parent.name,
                           it.parent.name.split('\\.')[1],
                           it.baseName.split('_')[0],
                           it.baseName.split('_')[1],
                           it.baseName.split('_')[2].split('-')[0],
                           it.baseName.split('_')[2].split('-')[1],
                           it.baseName.split('_')[0] + "_" + it.baseName.split('_')[2].split('-')[0] +"_"+ it.baseName.split('_')[1] +"_"+ it.baseName.split('_')[2].split('-')[1] +"_"+ "l" +  it.parent.name.split('\\.')[1],
                           it ) }.subscribe{ println it}
