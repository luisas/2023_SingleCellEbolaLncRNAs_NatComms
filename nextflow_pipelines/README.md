# Pipelines

All the pipelines used for the computations in the paper.
Scripts to run the pipelines in a LSF environment are in 00_CLUSTER_Nord3.
Scripts to run the pipelines in a SLURM environment are in 01_CLUSTER_Marenostrum.

Configurations files are in the configs folder.
Scripts used in the pipelines are in the scripts folder.

An example on how pipelines were executed by:
```
cd 01_CLUSTER_Marenostrum
sbatch 01_benchmak_annotation.batch
...etc..
```
