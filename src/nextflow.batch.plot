#!/bin/bash
#SBATCH --job-name=NextflowQAPipeline
#SBATCH --time=00-03:00:00
#SBATCH --workdir=.
#SBATCH --error=err/Nextflow-%j.err
#SBATCH --output=out/Nextflow-%j.out
#SBATCH --nodes=2
#SBATCH --cpus-per-task=48
#SBATCH --exclusive

module load nextflow
module load java
module load curl

# For conversion
module load java picard/2.20.0
# For quality assessment
module load fastqc/0.11.5
# For Mapping
module load gcc/7.2.0 hisat2 samtools

#For plotting the dag Graph
module load graphviz/2.40.1

srun nextflow run /gpfs/projects/bsc83/Ebola/code/ebola/src/quality_assesment.nf -with-dag flowcharts/flowchart_rnaseq_zyagen.png
