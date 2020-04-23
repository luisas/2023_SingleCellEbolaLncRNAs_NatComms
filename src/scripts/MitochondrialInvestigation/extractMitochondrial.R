shhh <- suppressPackageStartupMessages
shhh(library(rtracklayer))

args = commandArgs(trailingOnly=TRUE)
file = args[1]
#file = "/home/luisas/Desktop/cluster/gene_annotations/gencode/release_26/gencode.v26.basic.annotation.gtf"
output = args[2]

gtf<- import(file)
mt_gtf <- gtf[seqnames(gtf) == "MT"]

export(mt_gtf, file.path(output))
