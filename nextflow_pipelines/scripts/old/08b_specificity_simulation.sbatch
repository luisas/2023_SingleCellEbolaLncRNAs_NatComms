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


# In vivo
PREFIX="/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/04_simulation"
INPUT_SIMULATION=${PREFIX}"/m_seurat.rds"
DF_CELLTYPE=${PREFIX}"/df_celltype.rds"

#srun Rscript ${SCRIPT} ${INPUT_SIMULATION} ${PREFIX}"/01_specificity_fc/00_specificity_fc.rds" "-"
#srun Rscript ${SCRIPT_PROP} ${INPUT_SIMULATION} ${PREFIX}"/02_specificity_prop/01_specificity_prop.rds" "-"
#srun Rscript ${SCRIPT_CHI} ${DF_CELLTYPE} ${PREFIX}"/03_spec_chi/spec_chi.rds" ${INPUT_SIMULATION}
#srun Rscript ${SCRIPT_PROP_MOD} ${DF_CELLTYPE} ${PREFIX}"/04_spec_propmod/spec_propmod.rds" ${INPUT_SIMULATION}
#srun Rscript ${SCRIPT_TAU} ${DF_CELLTYPE} ${PREFIX}"/04_spec_propmod/spec_propmod.rds" ${INPUT_SIMULATION}
