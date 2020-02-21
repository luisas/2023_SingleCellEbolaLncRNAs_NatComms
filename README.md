# Ebola

Pipeline to obtain lncrna prediction was run in marenostrum with: 

sbatch /gpfs/projects/bsc83/Projects/Ebola/code/ebola/src/01_CLUSTER_Marenostrum/final_annotation_external.batch
sbatch /gpfs/projects/bsc83/Projects/Ebola/code/ebola/src/01_CLUSTER_Marenostrum/final_annotation_batchzyagen.batch

in nord3 with : 

bsub < /gpfs/projects/bsc83/Projects/Ebola/code/ebola/src/00_CLUSTER_Nord3/final_annotation.batch

Example in Marenostrum 
#!/bin/bash
#SBATCH --job-name=NextflowQAPipeline
#SBATCH --workdir=.
#SBATCH --error=err/Nextflow-%j.err
#SBATCH --output=out/Nextflow-%j.out
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48


module load java/8u131
module load intel/2017.1
module load nextflow/19.03.0

module load singularity/3.2.0


srun nextflow run /gpfs/projects/bsc83/Projects/Ebola/code/ebola/src/Annotation.nf \
                  --strandness "RF" \
                  --umi "false" \
                  --output_dir_name "02_RNA-Seq_external" \
                  --fastqs "/gpfs/projects/bsc83/Data/Ebola/00_RawData/extrenal_rhesus_RNA-Seq/*/*/*/*/*/*.{1,2}.fastq.gz" \
                  -w /gpfs/projects/bsc83/Data/Ebola/01_Ebola-RNASeq/work/ \
                  -c ../nextflow.config.rnaseq
