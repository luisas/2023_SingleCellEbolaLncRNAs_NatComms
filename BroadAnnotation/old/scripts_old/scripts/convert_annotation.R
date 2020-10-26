
#!/usr/bin/env Rscript
library(rtracklayer)
args = commandArgs(trailingOnly=TRUE)

file = args[1]
annotation = args[2]
out = args[3]


gtf <- import(file)
seqlevelsStyle(gtf) <- annotation
export(gtf,out)