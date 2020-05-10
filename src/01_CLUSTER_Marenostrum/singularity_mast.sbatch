#!/bin/bash
#SBATCH --job-name="singularity_mastDE"
#SBATCH --output="out/singularity_mastDE.%j.out"
#SBATCH --error="err/singularity_mastDE.%j.err"
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48

module load java/8u131
module load intel/2017.1

module load singularity/3.2.0

CONTAINER=/gpfs/projects/bsc83/utils/containers/singularity/rdeanalysis.simg

mpirun singularity exec ${CONTAINER} Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/src/scripts/SeuratForCluster/MAST_DE_ebolaremoved.R "immune.combined.rds"
