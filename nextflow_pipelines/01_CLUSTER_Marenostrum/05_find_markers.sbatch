#!/bin/bash
#SBATCH --job-name=NextflowQAPipeline
#SBATCH --workdir=.
#SBATCH --error=Nextflow-%j.err
#SBATCH --output=Nextflow-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --qos=debug


module load java/8u131
module load intel/2017.1
module load R/3.6.1

cd /gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/
chmod a+x 02_correlation_cellcyle.R
./02_correlation_cellcyle.R


#./02_identify.markers.R
