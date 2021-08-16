#!/bin/bash
#SBATCH --job-name=PreprocessZyagen
#SBATCH --workdir=.
#SBATCH --error=err/Nextflow-%j.err
#SBATCH --output=out/Nextflow-%j.out
#SBATCH --cpus-per-task=48
#SBATCH --qos=debug

module load bedtools 

bedtools getfasta -name -s -fi "/gpfs/projects/bsc83/Data/assemblies/ensembl/release-100/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa" -bed "/gpfs/projects/bsc83/Data/Ebola/03_tfbs/00_promoters.bed" -fo "/gpfs/projects/bsc83/Data/Ebola/03_tfbs/01_promoters.fa"
