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
# Boolean: should it execute the ranked corretion or not
correction <- args[3]


# Read in the Seurat object
immune.combined <- readRDS(file)


# ---------------------------------
#        Ranking Correction 
# ---------------------------------

if(correction=="corrected"){
  cell_correction_n_genes <-  function(column, threshold  = 1){
    column[column < threshold] <- NA
    ranking <- rank(as.vector(column), na.last = "keep")
    return(ranking/length(ranking))
  }
  corrected_matrix <- (apply(immune.combined@assays$RNA@data, 2,cell_correction_n_genes ))
  expression_matrix <- corrected_matrix
}else{
  expression_matrix <- immune.combined@assays$RNA@data
}

#############    ORDER THE MATRIX    ###################    
# Colnames: cell identities (cell-types)
colnames(expression_matrix) <- Idents(immune.combined)
# Rownames: genes
rownames(expression_matrix) <- rownames(immune.combined)
########################################################



# Calculate the mean 
weighted_mean <- function(expression_matrix, gene, ident, threshold=1){
  

  # Extract only the expression of the gene in the cell-type which we aer analyzing
  expression <- as.vector(expression_matrix[gene,colnames(expression_matrix) == ident])
  
  # Only retain cells in which the gene shows expression, assign NAs to the others
  expression[expression < threshold] <- NA

  
  # Calculate the mean expression
  #mean_expression <- ifelse(is.na(mean(expression, na.rm=T)),2, (mean(expression, na.rm=T)))
  mean_expression <- mean(expression, na.rm=T)
  
  # Calculate the Weight the mean expression of the gene by the proportion of cells in which it express 
  weight <-  sum(!is.na(expression))/length(expression)
  # Weight the mean:  Multiply the expression vector by the fraction of cells expressed in each cell tyoe
  return(mean_expression*weight)
  #return(mean_expression)
}


calc_score_gene <- function(expression_matrix, gene){
  # Calculate the mean expression of the gene in each cell-type (identities)
  identities <- unique(colnames(expression_matrix))
  weighted_means <- unlist(lapply(identities, function(ident) weighted_mean(expression_matrix,gene , ident )))
  
  # Create Fold-Changes Matrix
  res <- (outer(weighted_means, weighted_means, foldchange))
  diag(res) <- NA 
  colnames(res) <- rownames(res) <- identities
  
  # Per cell-type Retain the minimum out of the calculated fold-changes as the final score
  scores <- apply(res, 1, function(x) min(x, na.rm = T))
  # Use the maximum, acorss all celltype, as the final cell-type specificity score fot the gene 
  gene_score <- cbind(data.frame(gene=gene, score = max(as.vector(scores)), stringsAsFactors = F),t(scores))
  return(gene_score)
}



cl <- makeCluster(4, type="FORK")
clusterExport(cl,list("expression_matrix","weighted_mean","calc_score_gene"),envir=globalenv())

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
  clusterExport(cl,list("expression_matrix_permuted","weighted_mean","calc_score_gene"),envir=globalenv())
  
  spec_scores <- parLapply(cl, rownames(expression_matrix_permuted), function(gene) calc_score_gene(expression_matrix_permuted, gene))
  final <- Reduce("rbind", spec_scores)
  saveRDS(final, file.path(dirname(outfile), paste("Permutations_fc_",i,".rds", sep = "")))
  
  i <- i+1
}

stopCluster(cl)

