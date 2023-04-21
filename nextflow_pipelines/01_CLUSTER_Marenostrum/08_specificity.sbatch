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

SCRIPT_TAU=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity_tau.R
SCRIPT_PROP_MOD_RANGED=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/03_specificity_alternativescore.R



# IN VIVO
PREFIX="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects"
ININVIVO=${PREFIX}"/03_prep/03_immune.combined.ready.rds"
DF_CELLTYPE_INVIVO=${PREFIX}"/05_stats/df_celltype.rds"

srun Rscript ${SCRIPT_TAU} ${ININVIVO} ${PREFIX}"/05_stats/00_specificity/03_specificity_tau.rds"
srun Rscript ${SCRIPT_PROP_MOD_RANGED} ${DF_CELLTYPE_INVIVO} ${PREFIX}"/05_stats/00_specificity/04_specificity_alternativescore.rds"

# Replication: Prep
# 1. ---- Prep
# srun Rscript "/gpfs/projects/bsc83/Projects/Ebola/code/ebola/analysis/04_general_expression_patterns_lnc_vs_pc/01_prep33k.R"
# 2. ----
# SCRIPT=/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/06_calc_stats.R
# Rscript ${SCRIPT} "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/01_immunecombined_idents.rds" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/05_stats"
#
#
# PREFIX="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs"
# INREP=${PREFIX}"/01_immunecombined_idents.rds"
# DF_CELLTYPE_REP=${PREFIX}"/05_stats/df_celltype.rds"
#
# Rscript ${SCRIPT_TAU} ${INREP} ${PREFIX}"//03_specificity_tau.rds"
# Rscript ${SCRIPT_PROP_MOD_RANGED} ${DF_CELLTYPE_REP} ${PREFIX}"/04_specificity_alternativescore.rds"
