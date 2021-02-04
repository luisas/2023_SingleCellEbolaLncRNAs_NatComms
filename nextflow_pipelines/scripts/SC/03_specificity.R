library(Seurat)
library(gtools)


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
  if(correction!="corrected"){
    # Only retain cells in which the gene shows expression, assign NAs to the others
    expression[expression < threshold] <- NA
  }
  
  # Calculate the mean expression
  mean_expression <- mean(expression, na.rm=T)
  # Calculate the Weight the mean expression of the gene by the proportion of cells in which it express 
  weight <-  sum(!is.na(expression))/length(expression)
  # Weight the mean:  Multiply the expression vector by the fraction of cells expressed in each cell tyoe
  return(mean_expression*weight)
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



# Calculate the specificity score for all the genes
spec_scores <- lapply(rownames(expression_matrix), function(gene) calc_score_gene(expression_matrix, gene))
#spec_scores_min <- lapply(df_test_min$gene, function(gene) calc_score_gene(expression_matrix, gene))

# Save
final <- Reduce("rbind", spec_scores)
saveRDS(final, outfile)


