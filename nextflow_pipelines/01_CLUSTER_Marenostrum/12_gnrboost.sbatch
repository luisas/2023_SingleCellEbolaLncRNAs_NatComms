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


CONTAINER=/gpfs/projects/bsc83/utils/containers/singularity/aertslab-pyscenic-0.10.0.sif
SCRIPT="/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/grnboost2.py"


srun singularity exec ${CONTAINER} python3 ${SCRIPT} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/matrix_monocytes_de.csv" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/output_monocytes_de.tsv"
