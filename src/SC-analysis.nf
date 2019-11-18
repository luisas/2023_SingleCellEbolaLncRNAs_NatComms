// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"


params.dataset_bam_dir = "/gpfs/projects/bsc83/Data/Ebola/00_RawData/scRNAseq_exvivo_alin_bams"


dataset_bam = Channel
              .fromPath("${params.dataset_bam_dir}/*.bam")
              .ifEmpty('bam files directory is empty')
              .map { tuple(it.baseName.split('\\.')[0],
                           it.baseName.split('\\.')[1],
                           it.baseName.split('\\.')[3],
                           it.baseName.split('\\.')[6],
                           it.baseName.split('\\.')[0]+ "_" + it.baseName.split('\\.')[6] + "_" + it.baseName.split('\\.')[1] + "_" + it.baseName.split('\\.')[3],
                           it ) }.println{ it }



// process dedupUmi{
//   label 'big_mem'
//   cpus 24
//   storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"
//
//   input:
//   set animal_id, hpi, exp, replicate, complete_id, file(bam) from dataset_bam
//
//   output:
//   set dataset_name, dayPostInfection, tissue, sample, complete_id,
//       file("${complete_id}.UMI.f3.q60.umi_dedup.bam"),
//       file("${complete_id}.UMI.f3.q60.umi_dedup.bam.bai") into (dedup_umi_bams, dedup_umi_bams_for_count, dedup_umi_bams_for_distr, dedup_umi_bams_for_count_2, dedup_umis_3, dedup_umis_4)
//   //file "umi_logs" into umi_logs
//
//   script:
//   """
//   mkdir umi_logs
//   umi_tools dedup --stdin=${filtered_bam} --log=${complete_id}.umi_tools.dedup.log --output-stats=${complete_id} --extract-umi-method=tag --umi-tag=RX --paired --method directional  > ${complete_id}.UMI.f3.q60.umi_dedup.bam
//   mv ${complete_id}.umi_tools.dedup.log umi_logs/
//   mv ${complete_id}_per* umi_logs/
//   mv ${complete_id}_edit* umi_logs/
//   samtools index ${complete_id}.UMI.f3.q60.umi_dedup.bam
//   """
//
// }
