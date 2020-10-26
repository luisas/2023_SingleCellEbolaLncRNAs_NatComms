



log.info "=============================================="
log.info " Merge stringtie assemblies and perform prediction "
log.info "=============================================="


// BaseFolders
params.prefix_data = "/gpfs/projects/bsc83/Data"
params.output_dir = "${params.prefix_data}/Ebola/99_BroadAnnotation/"


// Reference annotation
params.gtf = "${params.prefix_data}/Ebola/99_BroadAnnotation/03_novel_lncRNAs_list/00_prefilter_candidates/prefilter_candidates.gtf"
our = Channel.fromPath("${params.gtf}")


files = Channel.fromPath("/gpfs/projects/bsc83/Data/Ebola/00_RawData/BroadTranscriptomesComplete/bams/*.bam")
                .ifEmpty("No bams found")
                .map { tuple(it.baseName, it) }


files.into{ bams; }



process stringtie{
   cpus 1
   storeDir "${params.output_dir}/04_requantification_broad_notstrand/"

   input:
   set file(gtf_names) from our.collect()
   set complete_id,
       file(bam) from bams

   output:
   file("*") into expression_channel

   script:
   """
   # Calc the counts for the umi_dedup
   stringtie -eB -G ${gtf_names} ${bam} -A ${complete_id}.gene_abundances.tsv
   mv t_data.ctab ${complete_id}_t_data.ctab
   """
}



 workflow.onComplete {
 	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
 }
