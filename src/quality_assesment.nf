#!/usr/bin/env nextflow

/*
 * Description
 */

params.dataset_bam_dir = "/gpfs/scratch/bsc83/bsc83024/test_dataset/bams_per_lane/*/"
params.output_dir = "/gpfs/scratch/bsc83/bsc83024/test_output/"

dataset_bam = Channel
              .fromPath("${params.dataset_bam_dir}/*_long.bam")
              .map { tuple(it.parent.name,
                           it.baseName.split('_')[1],
                           it.baseName.split('_')[1],
                           it.baseName.split('_')[2].split('-')[0],
                           it.baseName.split('_')[2].split('-')[1],
                           it ) }


process convert_bam_to_fastqc {

    publishDir "${params.output_dir}/01_fastq/$dataset_name/$tissue/$sample"

    input:
    set lane, dataset_name, dayPostInfection, tissue, sample, file(unmapped_bam) from dataset_bam

    output:
    set file("${sample}.1.fq"), file("${sample}.2.fq"), file('${sample}.unpaired.fq')

    script:
    """
    java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar SamToFastq I=${unmapped_bam} FASTQ=${sample}.1.fq  SECOND_END_FASTQ=${sample}.2.fq UNPAIRED_FASTQ=${sample}.unpaired.fq
    """

}
