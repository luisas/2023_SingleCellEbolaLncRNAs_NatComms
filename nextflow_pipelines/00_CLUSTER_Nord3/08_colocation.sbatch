#!/bin/bash
#BSUB -J NextflowannotationPipeline
#BSUB -cwd .
#BSUB -W 48:01
#BSUB -e err/Nextflow-%J.err
#BSUB -o out/Nextflow-%J.out

module load gcc/6.1.0
module load CURL/7.49.0
module load BZIP2/1.0.6

module load R/3.6.1


Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/08_colocation.R 

