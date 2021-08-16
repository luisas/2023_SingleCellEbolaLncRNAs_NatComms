

# ---------------- Prepare object for replication of expression patterns analysis of pc and lnc
# 1. create seurat object
# 2. Identity cell-types
# 3. Store informations for future analyses 

suppressPackageStartupMessages(library(rtracklayer))
suppressPackageStartupMessages(library(Seurat))

# Imports
clusterpath <- "/gpfs/projects/bsc83/Data/"
gene_annotation_path <- file.path(clusterpath, "gene_annotation")
data_path <- file.path(clusterpath, "Ebola")
plots <- file.path(data_path,"plots")
# Palette used throughout the script
col_lnc = "navy"
col_mrna = "#8DC3A7"
palette_plot_percentage <- c(col_lnc, col_mrna)
set.seed(123)
# Read reference files 
ref<- import(file.path(gene_annotation_path,"/ensembl_release100/homo_sapiens/Homo_sapiens.GRCh38.100.gtf"))
table(ref$gene_biotype)
lnc <- ref[ref$gene_biotype == "lncRNA",]$gene_name
length(unique(lnc))
pc <- ref[ref$gene_biotype == "protein_coding",]$gene_name

# Prepare Seurat object 
data <- Read10X(file.path(data_path,"02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/filtered_gene_bc_matrices/hg19/"))
immune.combined = CreateSeuratObject(counts = data)
dim(immune.combined)
immune.combined <- NormalizeData(immune.combined)
immune.combined <- FindVariableFeatures(immune.combined)
immune.combined <- ScaleData(immune.combined)
immune.combined <- RunPCA(immune.combined, npcs = 20, verbose = FALSE)
immune.combined <- RunUMAP(immune.combined, reduction = "pca", dims = 1:20)
immune.combined <- FindNeighbors(immune.combined, reduction = "pca",dims = 1:20 )
immune.combined <- FindClusters(immune.combined, resolution = 0.03)
saveRDS(immune.combined, file.path(data_path,"02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/01_immunecombined_clustered.rds"))


#immune.combined <- readRDS(file.path(data_path,"02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/01_immunecombined_clustered.rds"))

# Save information about gene biotypes
annotated_lncrnas <- rownames(immune.combined)[rownames(immune.combined) %in% lnc]
saveRDS(annotated_lncrnas, file.path(data_path, "02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/05_stats/all_lncrnas.rds"))
annotated_mrnas <- rownames(immune.combined)[rownames(immune.combined) %in% pc]
saveRDS(annotated_mrnas, file.path(data_path, "02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/05_stats/annotated_mrnas.rds"))


#DimPlot(immune.combined, reduction = "umap", label = TRUE, label.size = 9)+theme_void()

# remove ambigous cluster of < 150 cells
keep_cluster <- names(table(Idents(immune.combined))[as.vector(table(Idents(immune.combined)) > 150)])
immune.combined <- immune.combined[,Idents(immune.combined) %in% keep_cluster]
#DimPlot(immune.combined, reduction = "umap", label = TRUE, label.size = 9)+theme_void()

b <- c("CD79B", "MS4A1", "CD19", "IGHM")
CD8T <- c("CD3D", "GZMB", "GNLY")
CD4T <-  c("CD3D", "IL7R")
t <- c(CD4T, CD8T)
nk <- c("KLRB1", "GZMB","FCGR3")
mono <- c("LYZ", "PSAP", "CFD") 
neut <- c( "CD177","LCN2") 
marker.genes_red <- c(t,b,mono,neut)


#FeaturePlot(immune.combined, b, order = T )
#FeaturePlot(immune.combined, t, order = T )
#FeaturePlot(immune.combined, nk, order = T )
#FeaturePlot(immune.combined, mono, order = T )

immune.combined <- RenameIdents(immune.combined, `0` = "T", `1` = "Monocyte", `2` = "B", 
                                `3` = "NK")

#DimPlot(immune.combined, reduction = "umap", label = TRUE, label.size = 9)+theme_void()
saveRDS(immune.combined, file.path(data_path,"02_scRNA-Seq_PBMCs/03_validation_pbmcs_external/01_10xGenomic_10kPBMCs/01_immunecombined_idents.rds"))