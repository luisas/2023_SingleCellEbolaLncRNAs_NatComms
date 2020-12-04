
library(Seurat)
library(purrr)
library(stringr)
library(scater)
library(rtracklayer)
library(reshape2)
library(Matrix)

args = commandArgs(trailingOnly=TRUE)
#data_dir <- "/home/luisas/Desktop/cluster/data/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/04_DigitalExpressionMatrix_OK"
#outfile <- "PBMC_after_gene_QC.rds"
data_dir <- args[1]
ref <- import(args[2])
outfile_rds <- args[3]
outfile_mtx <- args[4]


iterate_files <- function(inpath, pattern_string){
  files <- list.files(path=inpath, pattern= pattern_string, full.names=TRUE, recursive=TRUE)
  return(files)
}


# ----------------------
# Obtain Seurat Objects 
# ----------------------
get_Seurat_object <- function(file){
  pbmc.data <- read.table(file = file, header = TRUE, row.names = 1)
  proj <-str_replace_all(str_replace_all(unlist(rev(unlist(str_split(file, pattern= "/"))[-1])[[1]]), "_", "-"), ".dge.txt.gz", "")
  proj <- str_sub(proj, 1, str_length(proj)-3)
  pbmc <- CreateSeuratObject(pbmc.data, project = proj, min.cells = 3, min.features = 200)
  return(pbmc)
}



# Load the PBMC dataset
files <- iterate_files(data_dir, ".dge.txt.gz")
length(unique(files_umis))
files_umis <- files[!str_detect(files, "reads" )]


##---------------------------------------------------------
#  1. MERGE quantification  
##---------------------------------------------------------
# all matrices files from single cell, still separated per lane
list <- map(files_umis[1:length(files_umis)], get_Seurat_object)

list_names <- lapply(files_umis, function(x) str_replace_all(unlist(rev(unlist(str_split(x, pattern= "/"))[-1])[[1]]), "_", "-"))
# Remove objects where nocells were identified 
list_detected_only <- list[unlist(lapply(list, function(x) dim(x)[1] != 0 ))]
list_names <- list_names[unlist(lapply(list, function(x) dim(x)[1] != 0 ))]

pbmc <- list_detected_only[1][[1]]
list_detected_only <- list_detected_only[-1]
# Merge all seurats object in one 
pbmc.big <- merge(pbmc, y = list_detected_only, add.cell.ids = list_names, project = "PBMCbig")
#saveRDS(pbmc.big, file.path(result_path,"EXVIVO_rhemac10_merged.rds"))

##---------------------------------------------------------
##                           Cell QC 
##  Remove cells with too little or too high UMIS count 
##  Remove cell based on number of detected genes per cell
##----------------------------------------------------------
#pbmc.big <- readRDS(file.path(result_path,"EXVIVO_rhemac10_merged.rds"))

pbmc_sc <- as.SingleCellExperiment(pbmc.big)
pbmc_sc <- calculateQCMetrics(pbmc_sc)

# Library size
#lib_size_hist <- plot_histogram( sce = pbmc_sc, qc_metric = "total_counts", title = "Library Size (total UMI)", log_scale = FALSE)
# Number of detected genes
#n_detec_hist <- plot_histogram(sce = pbmc_sc,  qc_metric = "total_features_by_counts",title = "Number of detected genes", log_scale = FALSE)
#lsizeplot <- lib_size_hist +
#  geom_vline(xintercept = c(1000, 10000), linetype = "dashed", color = "red")+ xlim(0,20000)
#lsizeplot <- lib_size_hist +
#  geom_vline(xintercept = c(1000, 10000), linetype = "dashed", color = "red") + xlim(0,2000)
#ndectplot <- n_detec_hist +
#  geom_vline(xintercept = c(600, 2000), linetype = "dashed", color = "red")

pbmc <- pbmc_sc
keep_cells <-
  pbmc$total_counts > 1000 & pbmc$total_counts < 10000 &
  pbmc$total_features_by_counts > 600 & pbmc$total_features_by_counts < 2000 

pbmc <- pbmc[, keep_cells]
#saveRDS(pbmc,file.path(result_path,"EXVIVO_rhemac10_aftercellQC.rds"))


##----------------------------
##           Gene QC 
##  Exclude genes not showing up in a least 10 cells
##----------------------------  

n_cells <- rowSums(as.matrix(counts(pbmc)) > 0)
pbmc <- pbmc[n_cells > 10, ]
# saveRDS(pbmc, file.path(result_path,"seurat_pbmc_rhemac10_merged_aftercellandgeneqc.rds")
#pbmc <- readRDS(file.path(result_path,"EXVIVO_rhemac10_aftercellQC_aftercellandgeneqc.rds"))
##----------------------------
##           Mitochondrial QC 
##----------------------------  
pbmc <- as.Seurat(pbmc)
mitochondrial_genes <- ref[seqnames(ref) == "MT",]$gene_name
mitochondrial_genes <- mitochondrial_genes[mitochondrial_genes %in% rownames(pbmc)]
if(length(mitochondrial_genes)>0){
  pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, features = mitochondrial_genes)
  pbmc <- subset(pbmc, percent.mt < 5)
}

saveRDS(pbmc, outfile_rds)
# Extract count matrix
m <-as.matrix(GetAssayData(object = pbmc, slot = "counts"))
sparse.gbm <- Matrix(m , sparse = T )
writeMM(obj = sparse.gbm, file=outfile_mtx)


