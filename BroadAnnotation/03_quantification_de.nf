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
params.output_dir = "${params.prefix_data}/Ebola/99_BroadAnnotation_Feb2021/"

// BAM files for expression
bams = Channel.fromPath("/gpfs/projects/bsc83/Data/Ebola/00_RawData/BroadTranscriptomesComplete/bams-new/*.bam")
                .ifEmpty("No bams found")
                .map { tuple(it.baseName, it) }

// bams = Channel.fromPath("/gpfs/projects/bsc83/Data/Ebola/00_RawData/BroadTranscriptomesComplete/bams-new/A0204_S146_L002_dedup.bam")
//                 .ifEmpty("No bams found")
//                 .map { tuple(it.baseName, it) }

// Gtf
gtf = Channel.fromPath("/gpfs/projects/bsc83/Data/Ebola/99_BroadAnnotation_Feb2021/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf").collect()


/*
*   SCRIPTS
*/
prepde_script = Channel.fromPath("${baseDir}/scripts/prepDE.py").collect()


//bams.view { "value: $it" }


// process htseq{
//    cpus 1
//    storeDir "${params.output_dir}/04_htseq/${complete_id}"
//
//    input:
//    file(gtf)
//    set complete_id,
//        file(bam) from bams
//
//    output:
//    file("*") into expression_channel
//
//    script:
//    """
//    htseq-count ${bam} ${gtf} > ${complete_id}.htseqcount.csv
//    """
// }


process expression_stringtie{

  storeDir "${params.output_dir}/06_quantification/${complete_id}"

  input:
  file(gtf)
  file(prepde_script)
  set complete_id,
      file(bam) from bams

  output:
  file("*") into quantification_channel

  script:
  """
  # Calc the counts for the umi_dedup
  stringtie -e --fr -G ${gtf} ${bam} -A ${complete_id}.gene_abundances.tsv -o ${complete_id}.gtf
  echo "${complete_id} ./${complete_id}.gtf" > samples.txt
  ./${prepde_script} -i samples.txt -g ${complete_id}_gene_counts.csv -t ${complete_id}_transcript_counts.csv
  """
}

 workflow.onComplete {
 	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
 }
