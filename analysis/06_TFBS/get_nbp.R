# Read command line
args = commandArgs(trailingOnly=TRUE)

fimo_updated_prom <- readRDS(args[1])
ouput <- (args[2])

library("GenomicRanges")
library(parallel)



create_grange <- function(gene_motifs){
  return(Reduce("c", (lapply(split(gene_motifs, gene_motifs$motif_id), function(i){ GRanges(seqnames = as.character(i$chr),
                                                                                            strand = i$strand,
                                                                                            ranges = IRanges(start = i$start_updated,
                                                                                                             end = i$stop_updated,
                                                                                                             names = i$motif_id))}))))
}

get_nb_motifs <- function(gene, fimo){ data.frame(gene = gene, nbp = sum(width(reduce(create_grange(fimo[fimo$gene_id == gene, ])))))}


cl <- makeCluster(detectCores())
clusterExport(cl,c('gene','fimo_updated_prom','get_nb_motifs', "create_grange"))
clusterEvalQ(cl, library("GenomicRanges"))
sumary_nbp <- parLapply(cl,unique(fimo$gene_id), function(gene) get_nb_motifs(gene, fimo_updated_prom))
stopCluster(cl)
summary_nbp <- Reduce("rbind",sumary_nbp )
saveRDS(summary_nbp,output)