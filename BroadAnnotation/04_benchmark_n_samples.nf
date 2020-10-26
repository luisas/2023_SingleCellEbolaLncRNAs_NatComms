#!/usr/bin/env nextflow

/*
*
* Nextflow pipeline for filtering and quantifying the novel lncRNAs.
*
*/
log.info "=============================================="
log.info " Quantification and Filtering  "
log.info "=============================================="

// ------------------------------------------------------------
// ------------ INPUT PARAMETERS ------------------------------
// ------------------------------------------------------------

params.prefix_data = "/gpfs/projects/bsc83/Data"
params.output_dir = "${params.prefix_data}/Ebola/99_Broad_benchmark/"
params.output_dir_preliminary = "${params.prefix_data}/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/"
// Create 2 different batch files
//params.htseqsense="yes"
//params.umis = "true"
//params.strandrule="1++,1--,2+-,2-+"
params.prefix_label = "ribodepleted"



//printing.subscribe{ println it }


params.data_dir = "${params.output_dir}/"
//params.lnc_novel_compared_to_ref= "${params.data_dir}/01_stringtie_assembly_merged/01_gffCompare/merged.annotated.gtf"
params.assembly ="${params.data_dir}/01_stringtie_assembly_merged/stringtie_merged_reference_guided.gtf"
params.reference_annotated = "${params.output_dir_preliminary}/gene_annotations/UCSC/rheMac10_EBOV-Kikwit_UCSC.gtf"
reference_assembly = Channel.fromPath("${params.output_dir_preliminary}/reference_assembly/rheMac10_EBOV-Kikwit_UCSC.fa")

// // SCRIPTS
extractNovel_script = Channel.fromPath("${baseDir}/scripts/03_concordant_and_merge.R").collect()
prefilter_script = Channel.fromPath("${baseDir}/scripts/01_prefilter_script.R").collect()
expression_script = Channel.fromPath("${baseDir}/scripts/02_filterExpressed.R").collect()
mrnas_predicted = Channel.fromPath("${params.mrnas_predicted}").collect()
//assembly = Channel.fromPath("${params.assembly}").collect()
gtf2bed = Channel.fromPath("${baseDir}/scripts/gtf2bed.R").collect()
gtf_ref_1 = Channel.fromPath("${params.reference_annotated}").collect()
ref_annot = Channel.fromPath("${params.reference_annotated}").collect()
ref_annot.into{ ref; ref1; ref2; }
cpatmodel= Channel.fromPath("${params.output_dir_preliminary}/Human_logitModel.RData").collect()
cpathexamer= Channel.fromPath("${params.output_dir_preliminary}/Human_Hexamer.tsv").collect()
cpat_files = Channel.fromPath("${params.output_dir_preliminary}/cpat*").collect()

gtf_ref_1.into{gtf_ref; gtf_ref_2; gtfChannel2}
// /*
// * Extract only the novel (unkown and antisense) lncRNAs compared to the reference.
// * Also remove non concordant predictions.
// */

params.output_dir_sub = "test"
output_dir_sub= "${params.output_dir_sub}"



params.inputdir = "/gpfs/projects/bsc83/Data/Ebola/00_RawData/BroadTranscriptomesComplete/stringtie2-gtfs/benchmark_subsets/tissue_effect/0*"

dir = Channel.fromPath( "${params.inputdir}", type: 'dir' )
              .ifEmpty("No dirs found")
              .map { tuple(it.parent.toString().split('/').reverse()[0], it.baseName,it) }

dir.into{dirs; pr}

pr.subscribe{ println "$it"}

process StringTie_Merge_Reference_Guided{

  cpus 1
  storeDir "${params.output_dir}/$benchmark_type/$benchmark_subdir"


  input:
  set benchmark_type, benchmark_subdir, file(stringtie_gtfs) from dirs
  file reference_gtf from ref2

  output:
  set benchmark_type, benchmark_subdir, file("stringtie_merged_reference_guided.gtf") into (merged_denovo_assmebly, merged_de_novo_assembly_2)

  script:
  """
  stringtie --merge --fr -p ${task.cpus} -o stringtie_merged_reference_guided.gtf -G ${reference_gtf} ${stringtie_gtfs}/*
  """
}


process gffCompare2{

  storeDir "${params.output_dir}//$benchmark_type/$benchmark_subdir/01_gffcompare/"

  input:
   set benchmark_type, benchmark_subdir,file(merged_gtf) from merged_denovo_assmebly
   file reference_gtf from ref1

   output:
   file("merged*") into gff_compare_output_channel2

   script:
   """
   gffcompare -R -r ${reference_gtf} -o merged ${merged_gtf}
   """
 }
