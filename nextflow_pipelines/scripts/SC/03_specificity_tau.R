library(Seurat)
library(gtools)
library(parallel)
library(tispec)

# -------------------------------------------------------------------------
# -------------------------------------------------------------------------
#      This script calculates the cell-type specificity score of a gene
# -------------------------------------------------------------------------
# -------------------------------------------------------------------------


# Read in argumants
args = commandArgs(trailingOnly=TRUE)
# Seurat object (pre-processed and normalized)
file <- args[1]
# Where the result file should be saved
outfile <- args[2]
dir.create(dirname(file.path(outfile)), showWarnings = FALSE)


# Read in the Seurat object
immune.combined <- readRDS(file)
expression_matrix <- immune.combined@assays$RNA@data

#############    ORDER THE MATRIX    ###################    
# Colnames: cell identities (cell-types)
colnames(expression_matrix) <- Idents(immune.combined)
# Rownames: genes
rownames(expression_matrix) <- rownames(immune.combined)
########################################################

mean_exp <- function(expression_matrix, gene, ident){
  # Extract only the expression of the gene in the cell-type which we aer analyzing
  expression <- as.vector(expression_matrix[gene,colnames(expression_matrix) == ident])
  # Calculate the mean expression. Consider zeros. 
  mean_expression <- mean(expression, na.rm=F)
  if(is.na(mean_expression)){ mean_expression <- 0}
  
  return(mean_expression)
}

get_gene_mean_expressions <- function(gene, m){
  idents <- unique(colnames(m))
  v <- data.frame(t(unlist(lapply(idents, function(ident) mean_exp(m, gene, ident)))))
  colnames(v) <- idents
  rownames(v) <- gene
  return(v)
}

# ----------------------------------------------
#    Calculate specificity score on cluster 
# -----------------------------------------------
cl <- makeCluster(4, type="FORK")
clusterExport(cl,list("expression_matrix","mean_exp","get_gene_mean_expressions"),envir=globalenv())

# Calculate mean expression per cell-type 
means <- parLapply(cl,rownames(expression_matrix), function(gene) get_gene_mean_expressions(gene, expression_matrix))
m_exp <-Reduce("rbind",means)
print("-------------------------")
print("   Means calculated      ")
print("-------------------------")

# Calcualate TAU
tauExp <- calcTau((m_exp)) 
saveRDS(tauExp, outfile)
stopCluster(cl)

