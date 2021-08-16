

# Replication 
# Calc stats
SCRIPT=/home/luisa/Desktop/cluster/projects/code/ebola/nextflow_pipelines/scripts/SC/06_calc_stats.R
Rscript ${SCRIPT} "/home/luisa/Desktop/cluster/data/02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/01_immunecombined_idents.rds" "/home/luisa/Desktop/cluster/data/02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/05_stats"


PREFIX="/home/luisa/Desktop/cluster/data//02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs"
INREP=${PREFIX}"/01_immunecombined_idents.rds"
DF_CELLTYPE_REP=${PREFIX}"/05_stats/df_celltype.rds"

#srun Rscript ${SCRIPT_TAU} ${INREP} ${PREFIX}"//03_specificity_tau.rds"
srun Rscript ${SCRIPT_PROP_MOD_RANGED} ${DF_CELLTYPE_REP} ${PREFIX}"/04_specificity_alternativescore.rds"


