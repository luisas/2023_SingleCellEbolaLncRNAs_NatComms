#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(SingleCellExperiment))
shhh(library(Seurat))
shhh(library(stringr))
shhh(library(dplyr))
options(future.globals.maxSize = 10000 * 1024^2)


# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#      This script summarizes basic gene metrics
#      such as:
#      - Number of cells in which the gene is expressed (above a threshold)
#      - Median expression value, across the cells in which shows expression 
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------




calc_mean_and_percentage_cell <- function(subset,ident = "", df, threshold = 1 ){
  
  # When not calculated on a specific cell-type, ident == ""
  # Otherwise it corresponds to the cell-type 
  # Only subsetting when doing the calculation on a specific cell-type
  if(ident != ""){
    subset <- subset(subset, idents = ident)
  }
  
  # Selecting normalized data 
  expression_matrix <- as.matrix(subset@assays$RNA@data)
  
  
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


# --------------------------------------
# Read the command line arguments 
# --------------------------------------
args = commandArgs(trailingOnly=TRUE)
file = args[1] # Seurat pre-processed object
robjectsdir <- args[2] # Directory where gene classification files are stored
immune.combined <- readRDS(file)
all_lncrnas <- readRDS(file.path(robjectsdir, "all_lncrnas.rds"))
annotated_mrnas <- readRDS(file.path(robjectsdir, "annotated_mrnas.rds"))


# ------------------------------------------------------
#  Calculate general metrics ( not cell-type specific)
# ------------------------------------------------------

# Calculate metrics for lncRNAS
df_lnc <- data.frame( gene_id = c(), meanexpr = c(), perc_cells_expressing = c(), maxexpr = c(), medianexpr = c(), var = c(), n_cells = c(), tot_cells = c())
invisible(calc_mean_and_percentage_cell(subset(immune.combined, features = all_lncrnas), "",df_lnc, 1))

# Calculate metrics for PC genes 
df_mrna<- data.frame( gene_id = c(), meanexpr = c(), perc_cells_expressing = c(), maxexpr = c(), medianexpr = c(), var = c(), n_cells = c(),tot_cells = c())
invisible(calc_mean_and_percentage_cell(subset(immune.combined, features = annotated_mrnas), "",df_mrna, 1))

# Save 
saveRDS(df_lnc, file.path(robjectsdir, "df_lnc.rds"))
saveRDS(df_mrna, file.path(robjectsdir, "df_mrna.rds"))


# ------------------------------------------------------
# Calculate the metrics in specific cell-types separetly.
# ------------------------------------------------------
# Calculate metrics for lncRNAS
df_lnc_celltype <- data.frame( gene_id = c(), meanexpr = c(), perc_cells_expressing = c(), maxexpr = c(), medianexpr = c(), var = c(), n_cells = c())
# Collect number of genes expressed per cell type
subset <- subset(immune.combined, features = all_lncrnas)
invisible(lapply(as.vector(unique(Idents(subset))), function(x) calc_mean_and_percentage_cell(subset, x, df_lnc_celltype)))

# Calculate metrics for PC genes 
df_mrna_celltype <- data.frame( gene_id = c(),meanexpr = c(), not_na_cells = c(), maxexpr = c(), medianexpr = c(), var = c(), var = c())
# Collect number of genes expressed per cell type
subset <- subset(immune.combined, features = annotated_mrnas)
invisible(lapply(as.vector(unique(Idents(subset))), function(x) calc_mean_and_percentage_cell(subset, x, df_mrna_celltype)))

df_lnc_celltype$type <- "lncRNA"
df_mrna_celltype$type <- "mRNA"
df_celltype <- rbind(df_lnc_celltype, df_mrna_celltype)

# Save
saveRDS(df_celltype, file.path(robjectsdir, "df_celltype.rds"))
saveRDS(df_lnc_celltype, file.path(robjectsdir, "df_celltype_lnc.rds"))
saveRDS(df_mrna_celltype, file.path(robjectsdir, "df_mrna_celltype.rds"))
