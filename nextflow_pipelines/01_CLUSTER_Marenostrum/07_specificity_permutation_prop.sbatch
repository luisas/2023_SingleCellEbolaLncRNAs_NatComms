#!/bin/bash
#SBATCH --job-name=NextflowQAPipeline
#SBATCH --workdir=.
#SBATCH --error=err/Nextflow-%j.err
#SBATCH --output=out/Nextflow-%j.out
#SBATCH --cpus-per-task=48
#SBATCH --qos=debug

module load java/8u131
module load intel/2017.1
module load R/3.6.1

module load singularity


srun Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity_prop.R "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready_noneut.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/00_specificity_permutation_prop/specificity_real.rds" "notcorrected"
