# Nextflow pipelines for de novo lncRNAs annotation and Single-cell preprocessing

All the pipelines were run on the Marenostrum cluster.
Specific input parameters can be found in the *.batch files in the 01_CLUSTER_Marenostrum folder.

Pipeline were executed by: 
```
cd 01_CLUSTER_Marenostrum 
sbatch 01_annotation.batch
sbatch 01b_filter.batch 
sbatch 02_annotation_external.batch

...etc..

```