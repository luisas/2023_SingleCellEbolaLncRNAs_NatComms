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

module load singularity/3.2.0

CONTAINER=/gpfs/projects/bsc83/utils/containers/singularity/rdeanalysis.simg

srun singularity exec ${CONTAINER} Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/07_correlation_viraload.R "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/myeloids.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload"

srun singularity exec ${CONTAINER} Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/07_correlation_viraload.R "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/myeloids_infected.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload_infected"

srun singularity exec ${CONTAINER} Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/07_correlation_viraload.R "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/myeloids_infected_24.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload_infected_24"

srun singularity exec ${CONTAINER} Rscript /gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/07b_correlation_viraload_2.R "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/myeloids_infected_24.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload_infected_24_window100"
