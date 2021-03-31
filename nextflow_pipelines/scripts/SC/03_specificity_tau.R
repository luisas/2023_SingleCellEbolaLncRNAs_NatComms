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
  # Calculate the mean expression
  mean_expression <- mean(expression, na.rm=F)
  if(is.na(mean_expression)){ mean_expression <- 0}
  # Weight the mean:  Multiply the expression vector by the fraction of cells expressed in each cell tyoe
  return(mean_expression)
  #return(mean_expression)
}


get_gene_mean_expressions <- function(gene, m){
  #df <- data.frame( value = unlist(lapply(identities, function(ident) mean_exp(m, gene, ident))), ident = identities, stringsAsFactors = F ) 
  return(unlist(lapply(unique(colnames(m)), function(ident) mean_exp(m, gene, ident))))
}





# ----------------------------------------------
#    Calculate specificity score on cluster 
# -----------------------------------------------
cl <- makeCluster(4, type="FORK")
clusterExport(cl,list("expression_matrix","mean_exp","get_gene_mean_expressions"),envir=globalenv())

# Real data
# Calculate mean expression per cell-type 
means <- parLapply(cl,rownames(expression_matrix), function(gene) get_gene_mean_expressions(gene, expression_matrix))
m_exp <-matrix(Reduce("cbind",means), ncol = length(unique(rownames(expression_matrix))), byrow = T)
colnames(m_exp) <- unique(colnames(expression_matrix))
rownames(m_exp) <- rownames(expression_matrix)
# Prepare data frame 
df <- as.data.frame(m_exp)
rownames(df) <- rownames(m)
tauExp <- calcTau(as.data.frame(m_exp)) 
saveRDS(tauExp, outfile)


# Permutations 
n_perm <- 2
i <- 1
while (i <= n_perm){
  expression_matrix_permuted <- expression_matrix
  colnames(expression_matrix_permuted) <- permute(colnames(expression_matrix_permuted))
  clusterExport(cl,list("expression_matrix","mean_exp","get_gene_mean_expressions"),envir=globalenv())
  
  # Real data
  # Calculate mean expression per cell-type 
  means <- parLapply(cl,rownames(expression_matrix_permuted), function(gene) get_gene_mean_expressions(gene, expression_matrix_permuted))
  m_exp <-matrix(Reduce("cbind",means), ncol = length(unique(rownames(expression_matrix_permuted))), byrow = T)
  colnames(m_exp) <- unique(colnames(expression_matrix_permuted))
  rownames(m_exp) <- rownames(expression_matrix_permuted)
  # Prepare data frame 
  df <- as.data.frame(m_exp)
  rownames(df) <- rownames(m)
  tauExp <- calcTau(as.data.frame(m_exp)) 
  saveRDS(tauExp, file.path(dirname(outfile), paste("Permutations_tau_",i,".rds", sep = "")))
  
  i <- i+1
}

stopCluster(cl)

