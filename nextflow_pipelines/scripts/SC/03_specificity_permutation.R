
library(Seurat)
library(gtools)
library(parallel)

args = commandArgs(trailingOnly=TRUE)
#file <- "/home/luisas/Desktop/cluster/data/RObjects_old/results/seurat_pbmc_rhemac10_merged_aftercellandgeneqc_afterScrublet_dim_red.rds"
file <- args[1]
outfile <- args[2]
correction <- args[3]


immune.combined <- readRDS(file)



if(correction=="corrected"){
  cell_correction_n_genes <-  function(column, threshold  = 1){
    column[column < threshold] <- NA
    ranking <- rank(as.vector(column), na.last = "keep")
    return(ranking/length(ranking))
  }
  corrected_matrix <- (apply(immune.combined@assays$RNA@data, 2,cell_correction_n_genes ))
  rownames(corrected_matrix) <- rownames(immune.combined)
  colnames(corrected_matrix) <- Idents(immune.combined)
  expression_matrix <- corrected_matrix
}else{
  expression_matrix <- immune.combined@assays$RNA@data
  colnames(expression_matrix) <- Idents(immune.combined)
  rownames(expression_matrix) <- rownames(immune.combined)
}


weighted_mean <- function(expression_matrix, gene, ident, threshold=1){
  expression <- as.vector(expression_matrix[gene,colnames(expression_matrix) == ident])
  if(correction!="corrected"){
    expression[expression < threshold] <- NA
  }
  mean_expression <- mean(expression, na.rm=T)
  # Calculate the Weight the mean expression of the gene by the proportion of cells in which it express 
  weight <-  sum(!is.na(expression))/length(expression)
  # Weight the mean:  Multiply the expression vector by the fraction of cells expressed in each cell tyoe
  return(mean_expression*weight)
}


calc_score_gene <- function(expression_matrix, gene){
  identities <- unique(colnames(expression_matrix))
  weighted_means <- unlist(lapply(identities, function(ident) weighted_mean(expression_matrix,gene , ident )))
  # Create Fold-Changes Matrix
  res <- (outer(weighted_means, weighted_means, foldchange))
  diag(res) <- NA 
  colnames(res) <- rownames(res) <- identities
  scores <- apply(res, 1, function(x) min(x, na.rm = T))
  gene_score <- cbind(data.frame(gene=gene, score = max(as.vector(scores)), stringsAsFactors = F),t(scores))
  return(gene_score)
}

cl <- makeCluster(4, type="FORK")
clusterExport(cl,list("expression_matrix","weighted_mean","calc_score_gene"),envir=globalenv())

# Real data
#spec_scores <- parLapply(cl, rownames(expression_matrix), function(gene) calc_score_gene(expression_matrix, gene))
#final <- Reduce("rbind", spec_scores)
#saveRDS(final, outfile)


# Permutations 
n_perm <- 100
i <- 1
while (i <= n_perm){
  expression_matrix_permuted <- expression_matrix
  colnames(expression_matrix_permuted) <- permute(colnames(expression_matrix_permuted))
  clusterExport(cl,list("expression_matrix_permuted","weighted_mean","calc_score_gene"),envir=globalenv())
  
  spec_scores <- parLapply(cl, rownames(expression_matrix_permuted), function(gene) calc_score_gene(expression_matrix_permuted, gene))
  final <- Reduce("rbind", spec_scores)
  saveRDS(final, file.path(dirname(outfile), paste("Permutations_",i,".rds", sep = "")))
  
  i <- i+1
}

stopCluster(cl)




