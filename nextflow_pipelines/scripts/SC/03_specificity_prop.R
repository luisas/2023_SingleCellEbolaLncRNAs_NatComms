
library(Seurat)
library(gtools)
library(parallel)

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




calc_score_gene <- function(expression_matrix, gene){
  threshold <- 1
  expression <- expression_matrix[gene,expression_matrix[gene, ]>1]
  score <- max(table(names(expression))/sum(length(expression)))
  # Use the maximum, acorss all celltype, as the final cell-type specificity score fot the gene 
  gene_score <- data.frame(gene=gene, score =score, stringsAsFactors = F)
  return(gene_score)
}




cl <- makeCluster(4, type="FORK")
clusterExport(cl,list("expression_matrix","calc_score_gene"),envir=globalenv())

# Real data
spec_scores <- parLapply(cl, rownames(expression_matrix), function(gene) calc_score_gene(expression_matrix, gene))
final <- Reduce("rbind", spec_scores)
saveRDS(final, outfile)


# Permutations 
n_perm <- 2
i <- 1
while (i <= n_perm){
  expression_matrix_permuted <- expression_matrix
  colnames(expression_matrix_permuted) <- permute(colnames(expression_matrix_permuted))
  clusterExport(cl,list("expression_matrix_permuted","calc_score_gene"),envir=globalenv())
  
  spec_scores <- parLapply(cl, rownames(expression_matrix_permuted), function(gene) calc_score_gene(expression_matrix_permuted, gene))
  final <- Reduce("rbind", spec_scores)
  saveRDS(final, file.path(dirname(outfile), paste("Permutations_prop_",i,".rds", sep = "")))
  
  i <- i+1
}

stopCluster(cl)


