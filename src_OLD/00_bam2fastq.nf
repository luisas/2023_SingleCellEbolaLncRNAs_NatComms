
log.info "=============================================="
log.info "          Convert unmapped bam to fastq.gz    "
log.info "=============================================="

params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/"

unmapped_bams = Channel
                  .fromPath("${params.dataset_bam_dir_zyagen}/*_long.bam")
                  .ifEmpty('bam files directory is empty')
                  .map{ tuple(it.baseName.split('_')[0] + "_" + it.baseName.split('_')[2].split('-')[0] +"_"+ it.baseName.split('_')[1] +"_"+ it.baseName.split('_')[2].split('-')[1] +"_"+ "l" +  it.parent.name.split('\\.')[1],
                              it)}

process bam2fastq{

  storeDir "${params.output_dir}/01_fastq/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample,
            file(summary),
            file(sam), file(bam),
            file(bam) from unmapped_and_mapped_bams

  output:
  set file("rheMac10_EBOV-Kikwit.fa.index"), file("feelnc_codpot_out") into coding_potentials

  script:
  """
  bedtools bamtofastq [OPTIONS] -i $bam -fq <FASTQ>
  """

}
