#!/bin/bash
#SBATCH --job-name=NextflowQAPipeline
#SBATCH --workdir=.
#SBATCH --error=err/Nextflow-%j.err
#SBATCH --output=out/Nextflow-%j.out
#SBATCH --cpus-per-task=48
module load java/8u131
module load intel/2017.1
module load R/3.6.1
module load singularity
#module spider ANACONDA/5.1.0_py3
#module load python
#source /gpfs/projects/bsc83/apps/python/bin/activate
CONTAINER=/gpfs/projects/bsc83/utils/containers/singularity/aertslab-pyscenic-0.10.0.sif

#srun singularity exec ${CONTAINER} python3 "/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/grnboost.py" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/matrix_final_myeloids_de_ebola.csv" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/output_myeloids_de_ebola.tsv"

#srun singularity exec ${CONTAINER} python3 "/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/grnboost.py" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/matrix_final_b_de_ebola.csv" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/output_b_de_ebola.tsv"

srun singularity exec ${CONTAINER} python3 "/gpfs/projects/bsc83/Projects/Ebola/code/ebola/nextflow_pipelines/scripts/SC/grnboost.py" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/matrix_final_t_de_ebola.csv" "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/output_t_de_ebola.tsv"
