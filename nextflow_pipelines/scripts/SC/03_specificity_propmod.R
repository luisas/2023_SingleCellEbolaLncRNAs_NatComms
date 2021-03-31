
library(Seurat)
library(gtools)
library(parallel)
shhh <- suppressPackageStartupMessages
#shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(dplyr))
shhh(library(MAST))

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#      This script calculates the cell-type specificity score of a gene
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


# Read in argumants
args = commandArgs(trailingOnly=TRUE)
# Seurat object (pre-processed and normalized)
df_celltype <- readRDS(args[1])
# Where the result file should be saved
outfile <- args[2]
file_seurat <- args[3]

dir.create(dirname(file.path(outfile)), showWarnings = FALSE)


calc_score_gene <- function(df_celltype, gene){
  df_celltype_gene <- df_celltype[df_celltype$gene_id == gene, ]
  df_celltype_gene$prop <- df_celltype_gene$tot_cells/sum(df_celltype_gene$tot_cells)
  # Calculate the expected 
  df_celltype_gene$expected <- sum(df_celltype_gene$n_cells)*df_celltype_gene$prop
  df_celltype_gene$obs_prop <- df_celltype_gene$observed/sum(df_celltype_gene$observed)
  calc_score_2  <- function(df_celltype_gene, ident){
    df_celltype_gene_celltype <- df_celltype_gene[df_celltype_gene$Celltype == ident,];  
    s <- (df_celltype_gene_celltype$observed -df_celltype_gene_celltype$expected)/df_celltype_gene_celltype$expected
    return(s)
  }
  
  calc_score_1  <- function(df_celltype_gene, ident){
    df_celltype_gene_celltype <- df_celltype_gene[df_celltype_gene$Celltype == ident,];  
    s <- (df_celltype_gene_celltype$obs_prop -df_celltype_gene_celltype$prop)
    return(s)
  }
  score <- max(unlist(lapply(identities, function(ident){calc_score_1(df_celltype_gene, ident)} ))) 
  df_score <- data.frame(gene= gene, score = score, stringsAsFactors = F)
  return(df_score)
}

cl <- makeCluster(4, type="FORK")
clusterExport(cl,list("df_celltype","calc_score_gene"),envir=globalenv())
# Real data
spec_scores <- parLapply(cl, as.character(unique((df_celltype$gene_id))), function(gene) calc_score_gene(df_celltype, gene))
final <- Reduce("rbind", spec_scores)
saveRDS(final, outfile)
stopCluster(cl)



# Permutation 
# First permute the celltype labels, then compute the stats, then submit the score

# TODO add input file seurat 

# Load expression matrix
immune.combined <- readRDS(file_seurat)

all_lncrnas <- readRDS(file.path(dirname(outfile), "all_lncrnas.rds"))
annotated_mrnas <- readRDS(file.path(dirname(outfile), "annotated_mrnas.rds"))

expression_matrix <- immune.combined@assays$RNA@data

#############    ORDER THE MATRIX    ###################    
# Colnames: cell identities (cell-types)
colnames(expression_matrix) <- Idents(immune.combined)
# Rownames: genes
rownames(expression_matrix) <- rownames(immune.combined)
########################################################



calc_mean_and_percentage_cell <- function(expression_matrix_permuted,ident = "", df, threshold = 1 ){
  
  # When not calculated on a specific cell-type, ident == ""
  # Otherwise it corresponds to the cell-type 
  # Only subsetting when doing the calculation on a specific cell-type
  if(ident != ""){
    expression_matrix <- expression_matrix_permuted[,colnames(expression_matrix_permuted)==ident]
  }
  
  # Selecting normalized data 
  # expression_matrix <- as.matrix(subset@assays$RNA@data)
  
  
  # -------------------- TOTAL NUMBER OF CELLS ----------------------
  # Cells are columns and genes are rows.
  # The total number of cells is the same for each gene. 
  # Tot_cells:= a vector as long as the number of genes present in the matrix. 
  # All values are the same: the total number of columns ( cells ) in the matrix. 
  tot_cells <- rep(ncol(expression_matrix), nrow(expression_matrix))
  
  
  # ---------------- NUMBER OF CELLS SHOWING EXPRESSION ----------------
  # Count the number of cells per gene showing some expression
  # Only select cells that are expressed (> threshold) 
  # Assign NAs to other cells, just for ease of calculation later
  # and compute the number of cells that show some expression 
  expression_matrix[expression_matrix < threshold] <- NA
  not_na_cells <- tot_cells- rowCounts(expression_matrix, value = NA)
  perc_cells_expressing <- not_na_cells/tot_cells
  
  
  # ---------------- CALCULATE MEAN, MEDIAN, MAX, MIN EXPRESSION and VARIANCE ------------
  # Per gene, calculate the above mentioned metrics, EXCLUSIVELY on the 
  # cells where the gene shows expression 
  meanexpr <- rowMeans(expression_matrix, na.rm = TRUE)
  medianexpr <- rowMedians(expression_matrix, na.rm = TRUE)
  maxexpr <- rowMaxs(expression_matrix, na.rm = TRUE)
  minexpr <- rowMins(expression_matrix, na.rm = TRUE)
  var <- rowVars(expression_matrix, na.rm = TRUE)
  
  # Vector of all gene ids 
  gene_id <- rownames(expression_matrix)
  
  # Save 
  dfo <- deparse(substitute(df))
  df_new <- data.frame(gene_id = gene_id, meanexpr = meanexpr,medianexpr = medianexpr, perc_cells_expressing = perc_cells_expressing, maxexpr = maxexpr, var = var, n_cells = not_na_cells, tot_cells
                       = tot_cells)
  df_new$minexpr <- minexpr
  df_new$ident <- as.factor(ident)
  df_complete <- rbind(df, df_new)
  assign(dfo, df_complete, envir = .GlobalEnv);
}

cl <- makeCluster(4, type="FORK")

# Permutations 
n_perm <- 1
i <- 1
while (i <= n_perm){
  expression_matrix_permuted <- expression_matrix
  colnames(expression_matrix_permuted) <- permute(colnames(expression_matrix_permuted))
  expression_matrix_permuted <- as.matrix(expression_matrix_permuted)

  # Calculate metrics for lncRNAS
  df_lnc_celltype <- data.frame( gene_id = c(), meanexpr = c(), perc_cells_expressing = c(), maxexpr = c(), medianexpr = c(), var = c(), n_cells = c())
  # Collect number of genes expressed per cell type
  expression_matrix_permuted_lnc <- expression_matrix_permuted[all_lncrnas, ]
  invisible(lapply(as.vector(unique(colnames(expression_matrix_permuted_lnc))), function(x) calc_mean_and_percentage_cell(expression_matrix_permuted_lnc, x, df_lnc_celltype)))
  
  # Calculate metrics for PC genes 
  df_mrna_celltype <- data.frame( gene_id = c(),meanexpr = c(), not_na_cells = c(), maxexpr = c(), medianexpr = c(), var = c(), var = c())
  # Collect number of genes expressed per cell type
  expression_matrix_permuted_pc <- expression_matrix_permuted[annotated_mrnas, ]
  invisible(lapply(as.vector(unique(colnames(expression_matrix_permuted_pc))), function(x) calc_mean_and_percentage_cell(expression_matrix_permuted_pc, x, df_mrna_celltype)))
  
  # Save stats of permtuation
  saveRDS(df_lnc_celltype, file.path(dirname(outfile), "df_celltype_lnc_permutation.rds"))
  saveRDS(df_mrna_celltype, file.path(dirname(outfile), "df_mrna_celltype_permutation.rds"))
  df_celltype <- rbind(df_mrna_celltype,df_lnc_celltype )
  saveRDS(df_celltype, file.path(dirname(outfile),"df_celltype_permutation.rds"))
  
  
  # Export what is needed on cluster 
  clusterExport(cl,list("expression_matrix_permuted","calc_score_gene", "df_celltype"),envir=globalenv())
  
  
  
  spec_scores <- parLapply(cl, unique(df_celltype$gene_id), function(gene) calc_score_gene(df_celltype, gene))
  final <- Reduce("rbind", spec_scores)
  saveRDS(final, file.path(dirname(outfile), paste("Permutations_chi_",i,".rds", sep = "")))
  
  i <- i+1
}


stopCluster(cl)


