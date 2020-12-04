#!/usr/bin/env Rscript
shhh <- suppressPackageStartupMessages
shhh(library(Seurat))
shhh(library(stringr))
options(future.globals.maxSize = 10000 * 1024^2)


args = commandArgs(trailingOnly=TRUE)

file = args[1]
robjectsdir <- "/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjectsOLD/06_correlation/01_pearson"
immune.combined <- readRDS(file)

# Only calculate for the DE genes 
de_all_genes<- readRDS("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjectsOLD/04_DE/de_all_genes.rds")
immune.combined <- subset(immune.combined, features=de_all_genes)



identities <- unique(Idents(immune.combined))






get_correlation_matrix <- function(immune.combined, celltype, all = F){

  # Obtain the subset from the seura object and extract count matrix
  if(all == T){
    subset <- immune.combined
    celltype <- "all"
  }else{
    subset <- subset(immune.combined, idents = celltype)
  }

  count_matrix <- subset@assays$RNA@counts
  
  genes <- rownames(count_matrix)
  # Calculate rho and p-value from correlation 
  i = 1 
  j = 1 
  m_pval <- matrix(NA, length(genes), length(genes))
  m_rho<- matrix(NA, length(genes), length(genes))
  while (i <= length(genes)){
    j <- i+1
    while (j <= length(genes)){
      cor_object <- calc_correlation(genes[i],genes[j],  count_matrix, type = "pearson")
      if(is.na(cor_object)){
        m_pval[i,j] <- m_pval[j,i] <- NA
        m_rho[i,j] <- m_rho[j,i] <- NA
      }else{
        m_pval[i,j] <- m_pval[j,i] <- cor_object$p.value
        m_rho[i,j] <- m_rho[j,i] <- cor_object$estimate
      }
      j <- j+1
    }
    print(paste("Calculated gene: ", (genes[i]), sep = ""))
    i <- i+1
  }
  # Set colnames and rownames of the matrix
  rownames(m_pval) <- colnames(m_pval) <- rownames(m_rho) <- colnames(m_rho) <- genes
  
  saveRDS(m_pval, file.path(robjectsdir, paste("m_pval_",celltype,".rds", sep="")))
  saveRDS(m_rho, file.path(robjectsdir, paste("m_rho",celltype,".rds", sep="")))
}



identities <- Idents(immune.combined)
countmatrix <- immune.combined@assays$RNA@counts
gene1 <- "AAAS"
gene2 <- "AAAS"
gene1_vector <- as.vector(countmatrix[gene1,])
gene2_vector <- as.vector(countmatrix[gene2,])




calc_correlation_vectors <- function(gene1_vector,gene2_vector, type = "pearson"){
  mask <- gene1_vector !=0 & gene2_vector !=0
  # calculate the correlation coefficient
  if(sum(mask)>50){ 
    #pc <- bayes.cor.test(c1, c2) 
    pc <- cor.test(gene1_vector[mask], gene2_vector[mask], method = c(type))
    df <- data.frame( pval=pc$p.value, rho=pc$estimate)
    return(df)
  }else{ 
    return(NA)
  }
}

get_correlation_df_genes <- function(gene1_vector, gene2_vector, gene1, gene2, identities, identity){
   mask <- identities == identity 
   df <- calc_correlation(gene1_vector[mask], gene2_vector[mask])
   if(!is.na(df)){
     df$celltype <- identity
     df$g1 <- gene1
     df$g2 <- gene2
     return(df)
   }
}


cl <- makeCluster(4, type="FORK")


do.call(rbind,lapply(as.character(unique(identities)), function(x) get_correlation_df_genes(gene1_vector, gene2_vector, gene1, gene2, identities, x)))




expand.grid(rownames(countmatrix), rownames(countmatrix))

correlation <- parLapply(cl, rownames(expression_matrix_permuted), function(gene) calc_score_gene(expression_matrix_permuted, gene))


cl <- makeCluster(4, type="FORK")
clusterExport(cl,list("expression_matrix","weighted_mean","calc_score_gene"),envir=globalenv())



stopCluster(cl)









