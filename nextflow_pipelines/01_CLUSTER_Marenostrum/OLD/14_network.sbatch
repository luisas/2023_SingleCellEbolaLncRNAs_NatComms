#!/bin/bash
#SBATCH --job-name=NextflowQAPipeline
#SBATCH --workdir=.
#SBATCH --error=err/Nextflow-%j.err
#SBATCH --output=out/Nextflow-%j.out
#SBATCH --cpus-per-task=48
module load java/8u131
module load intel/2017.1
module load R/3.6.1
module load singularity/3.2.0
#module spider ANACONDA/5.1.0_py3
CONTAINER=/gpfs/projects/bsc83/utils/containers/singularity/rde.simg
#conda activate /gpfs/projects/bsc83/Data/Ebola/envs/arboreto-env

#srun Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/src/scripts/SeuratForCluster/geneQC.R "/gpfs/projects/bsc83/Projects/Ebola/code/ebola/src/scripts/Analysis/results/seurat_pbmc_rhemac10_merged_aftercellqc.rds"
#python3 "/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/grnboost.py"
srun singularity exec ${CONTAINER} Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/09_genie3.R
