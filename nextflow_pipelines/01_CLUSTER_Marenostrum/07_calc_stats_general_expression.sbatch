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

CONTAINER=/gpfs/projects/bsc83/utils/containers/singularity/rdeanalysis.simg
SCRIPT=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/06_calc_stats.R
SCRIPT_WITHZEROS=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/06_calc_stats_withzeros.R

srun singularity exec ${CONTAINER} Rscript ${SCRIPT} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats"
srun singularity exec ${CONTAINER} Rscript ${SCRIPT} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/05_stats"
#srun singularity exec ${CONTAINER} Rscript ${SCRIPT} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/v2_ut.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/05_stats"

srun singularity exec ${CONTAINER} Rscript ${SCRIPT_WITHZEROS} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats"
srun singularity exec ${CONTAINER} Rscript ${SCRIPT_WITHZEROS} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/05_stats"
