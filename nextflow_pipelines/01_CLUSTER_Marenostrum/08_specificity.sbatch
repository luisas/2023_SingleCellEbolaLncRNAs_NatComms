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

SCRIPT=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity.R
SCRIPT_PROP=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity_prop.R
SCRIPT_CHI=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity_chi_simulations.R
SCRIPT_PROP_MOD=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity_propmod.R
SCRIPT_TAU=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity_tau.R
SCRIPT_PROP_MOD_RANGED=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity_propmod_ranged.R


#srun Rscript ${SCRIPT_CHI} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/df_celltype.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/spec_chi.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/03_prep/03_immune.combined.ready.rds"

# IN VIVO
PREFIX="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects"
ININVIVO=${PREFIX}"/03_prep/03_immune.combined.ready.rds"
DF_CELLTYPE_INVIVO=${PREFIX}"/05_stats/df_celltype.rds"
DF_CELLTYPE_INVIVO_PERM=${PREFIX}"/05_stats/df_celltype_permutation.rds"
#srun Rscript ${SCRIPT_TAU} ${ININVIVO} ${PREFIX}"/05_stats/00_specificity/03_specificity_tau.rds"
srun Rscript ${SCRIPT_PROP_MOD_RANGED} ${DF_CELLTYPE_INVIVO} ${PREFIX}"/05_stats/00_specificity/04_specificity_propmod_ranged.rds" ${DF_CELLTYPE_INVIVO_PERM}

#srun Rscript ${SCRIPT} ${ININVIVO} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/00_specificity/00_specificity_fc.rds" "notcorrected"
#srun Rscript ${SCRIPT_PROP} ${ININVIVO} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/00_specificity/00_specificity_prop.rds"


# RUN it ex vivo
PREFIX="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects"
INEXVIVO=${PREFIX}"/03_prep/03_immune.combined.ready.rds"
srun Rscript ${SCRIPT_TAU} ${INEXVIVO} ${PREFIX}"/05_stats/00_specificity/03_specificity_tau.rds"
#srun Rscript ${SCRIPT} ${INEXVIVO} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/05_stats/00_specificity_fc/00_specificity_fc.rds" "notcorrectedS"
#srun Rscript ${SCRIPT_PROP} ${INEXVIVO} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/05_stats/00_specificity/01_specificity_prop.rds" "notcorrectedS"


# RUN CHI
#srun Rscript ${SCRIPT_CHI} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/df_celltype.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/spec_chi.rds" ${ININVIVO}
#srun Rscript ${SCRIPT_CHI} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/df_celltype_permutation.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/05_stats/spec_chi_permutation.rds"
#srun Rscript ${SCRIPT_CHI} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/05_stats/df_celltype.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/00_scRNA-Seq_exVivo_rhemac10/05_RObjects/05_stats/00_specificity/00_chi/spec_chi.rds" ${INEXVIVO}
