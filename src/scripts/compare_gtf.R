library(rtracklayer)
library(VennDiagram)

dir<- "/home/luisas/Desktop/cluster/proj/data"
lncRNAs_ref <- import(file.path(dir, "01_PreliminaryFiles/gene_annotations/rheMac8_EBOV-Kikwit_lncrna.gtf"))
lncRNAs_out <- import(file.path(dir, "02_RNA-Seq/05_lncrnaAnnotation/feelnc/feelnc_codpot_out/candidate_lncRNA.gtf.lncRNA.gtf"))

length(lncRNAs_out$transcript_id)
length(lncRNAs_ref$transcript_id)

library(RColorBrewer)
myCol <- brewer.pal(2, "Pastel2")
# Chart
vd_draw(list(lncRNAs_ref$gene_id, lncRNAs_out$gene_id))
vd_draw(list(lncRNAs_ref$transcript_id, lncRNAs_out$transcript_id))


lncRNAs_ref$transcript_id
lncRNAs_out$transcript_id

# ----- UTILS ---------------
vd_draw <- function(list){
  temp <- venn.diagram(
    x = list,
    category.names = c("Ref" , "FeelNC" ),
    filename = NULL,
    output=TRUE,
    
    # Circles
    lwd = 2,
    lty = 'blank',
    fill = c("#B3E2CD","#FDCDAC"),
    
    # Numbers
    cex = .7,
    fontface = "bold",
    fontfamily = "sans"
    
  )
  grid.draw(temp)
}


