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

SCRIPT=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/07_correlation_viraload.R
SCRIPT_PREP=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/00_removeEbolaReads.R

OUTPUT_INVIVO="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/06_correlation/03_viralload_infected_late"
OUTPUT_EXVIVO="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/06_correlation/03_viralload_infected_late"

PREFIX_INVIVO="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects"
PREFIX_EXVIVO="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects"



# ------------------------------
#         EX VIVO
# ------------------------------

# --- 1. Remove Ebola reads --------
 srun singularity exec ${CONTAINER} Rscript ${SCRIPT_PREP} ${PREFIX_EXVIVO}"/03_prep/03_immune.combined.ready.rds" \
                                                             ${PREFIX_EXVIVO}"/03_prep/ebolagenes.rds" \
                                                             ${PREFIX_EXVIVO}"/03_prep/03_immune.combined.ready_infectionstatus.rds" \
                                                            ${PREFIX_EXVIVO}"/03_prep/df_viralpercentage.rds" \
                                                            ${PREFIX_EXVIVO}"/03_prep/03_immune.combined.ready_infectionstatus_noebola.rds"

# Run correlations without EBOV reads
 srun singularity exec ${CONTAINER} Rscript ${SCRIPT} ${PREFIX_EXVIVO}"/03_prep/03_immune.combined.ready_infectionstatus_noebola.rds" \
                                                      "exvivo" \
                                                      ${OUTPUT_EXVIVO}"/mono_infected_late_ebolaremoved_correlation.rds"\
                                                      ${OUTPUT_EXVIVO}"/mono_infected_late_ebolaremoved.rds"





# # ------------------------------
#         IN VIVO
# ------------------------------
#srun singularity exec ${CONTAINER} Rscript ${SCRIPT_PREP} ${PREFIX_INVIVO}"/03_prep/03_immune.combined.ready.rds" \
#                                                           ${PREFIX_INVIVO}"/03_prep/ebolagenes.rds" \
#                                                           ${PREFIX_INVIVO}"/03_prep/03_immune.combined.ready_infectionstatus.rds" \
#                                                           ${PREFIX_INVIVO}"/03_prep/df_viralpercentage.rds" \
#                                                           ${PREFIX_INVIVO}"/03_prep/03_immune.combined.ready_infectionstatus_noebola.rds"


#srun singularity exec ${CONTAINER} Rscript ${SCRIPT} ${PREFIX_INVIVO}"/03_prep/03_immune.combined.ready_infectionstatus_noebola.rds" \
#                                                     "invivo" \
#                                                     ${OUTPUT_INVIVO}"/mono_infected_late_ebolaremoved_correlation.rds" \
#                                                     ${OUTPUT_INVIVO}"/mono_infected_late_ebolaremoved.rds"

