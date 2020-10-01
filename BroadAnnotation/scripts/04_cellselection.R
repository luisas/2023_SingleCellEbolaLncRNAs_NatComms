
library(readr)
library(dropbead)

# Read command line
args = commandArgs(trailingOnly=TRUE)

file <- args[1]
outfile <- args[2]
df_cell_readcount=read.table(file, header=F, stringsAsFactors=F)
# Estimate the number of cell
# threshold as in : https://doi.org/10.1016/j.celrep.2018.11.003
max.cells <- min(25000, nrow(df_cell_readcount))

#pdf("cumulative_distribution.pdf")
#p1 <- plotCumulativeFractionOfReads(df_cell_readcount, cutoff = max.cells)
#dev.off()
#print(c("Number or STAMPS:", estimateCellNumber(df_cell_readcount[1:max.cells, 1], max.cells)))

number_cells <- estimateCellNumber(df_cell_readcount[1:max.cells, 1], max.cells)
barcodes <- head(df_cell_readcount[order(-df_cell_readcount$V1),], n = number_cells)$V2

out<-file(outfile)
writeLines(barcodes, out)
close(out)
