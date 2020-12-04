#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(dplyr))
shhh(library(MAST))
options(future.globals.maxSize = 10000 * 1024^2)



calc_mean_and_percentage_cell <- function(subset,ident = "", df, threshold = 1 ){
  if(ident != ""){
    subset <- subset(subset, idents = ident)
  }
  # Using normalized data 
  expression_matrix <- as.matrix(subset@assays$RNA@data)
  
  # The totlal number of cells is the same for each gene. Cells are columns and genes are rows.
  tot_cells <- rep(ncol(expression_matrix), nrow(expression_matrix))
  
  # Count the number of cells per gene showing some expression
  # Only select cells that are expressed (> threshold) and say they are NA
  # and compute the number of cells that show some expression 
  expression_matrix[expression_matrix < threshold] <- NA
  meanexpr <- rowMeans(expression_matrix, na.rm = TRUE)
  medianexpr <- rowMedians(expression_matrix, na.rm = TRUE)
  maxexpr <- rowMaxs(expression_matrix, na.rm = TRUE)
  minexpr <- rowMins(expression_matrix, na.rm = TRUE)
  var <- rowVars(expression_matrix, na.rm = TRUE)
  not_na_cells <- tot_cells- rowCounts(expression_matrix, value = NA)
  perc_cells_expressing <- not_na_cells/tot_cells
  gene_id <- rownames(subset)
  
  # Save 
  dfo <- deparse(substitute(df))
  df_new <- data.frame(gene_id = gene_id, meanexpr = meanexpr,medianexpr = medianexpr, perc_cells_expressing = perc_cells_expressing, maxexpr = maxexpr, var = var, n_cells = not_na_cells, tot_cells
                       = tot_cells)
  df_new$minexpr <- minexpr
  df_new$ident <- as.factor(ident)
  df_complete <- rbind(df, df_new)
  assign(dfo, df_complete, envir = .GlobalEnv);
}



args = commandArgs(trailingOnly=TRUE)

file = args[1]
robjectsdir <- args[2]
immune.combined <- readRDS(file)

all_lncrnas <- readRDS(file.path(robjectsdir, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(file.path(robjectsdir, "annotated_mrnas.rds"))

df_lnc <- data.frame( gene_id = c(), meanexpr = c(), perc_cells_expressing = c(), maxexpr = c(), medianexpr = c(), var = c(), n_cells = c(), tot_cells = c())
invisible(calc_mean_and_percentage_cell(subset(immune.combined, features = all_lncrnas), "",df_lnc, 1))
df_mrna<- data.frame( gene_id = c(), meanexpr = c(), perc_cells_expressing = c(), maxexpr = c(), medianexpr = c(), var = c(), n_cells = c(),tot_cells = c())
invisible(calc_mean_and_percentage_cell(subset(immune.combined, features = annotated_mrnas), "",df_mrna, 1))

saveRDS(df_lnc, file.path(robjectsdir, "df_lnc.rds"))
saveRDS(df_mrna, file.path(robjectsdir, "df_mrna.rds"))


df_lnc_celltype <- data.frame( gene_id = c(), meanexpr = c(), perc_cells_expressing = c(), maxexpr = c(), medianexpr = c(), var = c(), n_cells = c())
# Collect number of genes expressed per cell type
subset <- subset(immune.combined, features = all_lncrnas)
invisible(lapply(as.vector(unique(Idents(subset))), function(x) calc_mean_and_percentage_cell(subset, x, df_lnc_celltype)))

# mRNA
df_mrna_celltype <- data.frame( gene_id = c(),meanexpr = c(), not_na_cells = c(), maxexpr = c(), medianexpr = c(), var = c(), var = c())
# Collect number of genes expressed per cell type
subset <- subset(immune.combined, features = annotated_mrnas)
invisible(lapply(as.vector(unique(Idents(subset))), function(x) calc_mean_and_percentage_cell(subset, x, df_mrna_celltype)))

saveRDS(df_lnc_celltype, file.path(robjectsdir, "df_celltype_lnc.rds"))
saveRDS(df_mrna_celltype, file.path(robjectsdir, "df_mrna_celltype.rds"))
