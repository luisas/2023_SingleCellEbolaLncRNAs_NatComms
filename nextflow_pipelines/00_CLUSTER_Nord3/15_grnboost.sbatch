#!/bin/bash
#BSUB -J NextflowannotationPipeline
#BSUB -cwd .
#BSUB -W 48:01
#BSUB -e err/Nextflow-%J.err
#BSUB -o out/Nextflow-%J.out



module load java/1.8.0u66
module load intel/2017.1
module load nextflow/19.03.0

module load singularity/3.2.0
#module spider ANACONDA/5.1.0_py3

CONTAINER=/gpfs/projects/bsc83/utils/containers/singularity/aertslab-pyscenic-0.10.0.sif

singularity exec ${CONTAINER} python3 "/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/grnboost.py"
