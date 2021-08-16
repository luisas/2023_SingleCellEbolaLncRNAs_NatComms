# Read command line
args = commandArgs(trailingOnly=TRUE)

fimo <- readRDS(args[1])
output <- args[2]

library("GenomicRanges")
library(parallel)
library(dplyr)


create_grange <- function(gene_motifs){
  return(Reduce("c", (lapply(split(gene_motifs, gene_motifs$motif_id), function(i){ GRanges(seqnames = as.character(i$chr),
                                                                                            strand = i$strand,
                                                                                            ranges = IRanges(start = i$start_updated,
                                                                                                             end = i$stop_updated,
                                                                                                             names = i$motif_id))}))))
}

get_max_cluster <- function(gene, fimo){
  r <- create_grange(fimo[fimo$gene_id == gene, ])
  df<- as.data.frame(findOverlaps(r))
  # Identify all clusters ( reduntant but we don't care )
  # Which is the size of the cluster to which that motif belongs to
  clusters <- df %>% dplyr::group_by(queryHits) %>%  dplyr::summarise(n = n())
  colnames(clusters) <- c("queryHits", "sizeCluster")
  return(data.frame(gene_id= gene, max_cluster =max(clusters$sizeCluster)))
}

library(parallel)
cl <- makeCluster(detectCores())
clusterExport(cl,c('fimo','get_max_cluster', "create_grange"))
clusterEvalQ(cl, library("GenomicRanges"))
clusterEvalQ(cl, library("dplyr"))
sumary_mc <- parLapply(cl,unique(fimo$gene_id), function(gene) get_max_cluster(gene, fimo))
stopCluster(cl)
summary_mc <- Reduce("rbind",sumary_mc )
saveRDS(summary_mc,output)