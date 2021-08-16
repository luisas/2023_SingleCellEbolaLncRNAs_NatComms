#!/bin/bash
#SBATCH --job-name=PreprocessZyagen
#SBATCH --workdir=.
#SBATCH --error=err/Nextflow-%j.err
#SBATCH --output=out/Nextflow-%j.out
#SBATCH --cpus-per-task=48
#SBATCH --qos=debug

module load R 

Rscript get_nbp.R "/gpfs/projects/bsc83/Data/Ebola/03_tfbs/fimo.rds" "/gpfs/projects/bsc83/Data/Ebola/03_tfbs/summary_nbp.rds"
