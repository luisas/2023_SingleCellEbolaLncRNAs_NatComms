
library(rtracklayer)
library(dplyr)

# Read command line
args = commandArgs(trailingOnly=TRUE)
cpc2 <- read.table(args[1])
cpat <- read.table(args[2], header = T ) 
cnit <- read.table(args[3], header = TRUE,sep="\t")
candidates <- import(args[4])
ref_gtf <- import(args[5])
outfile <- args[6]
outfile_name <- args[7]

#ribodepl_novel_compared_to_polya <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/03_polyA_vs_ribodepl/ribodeplVSpolyA.annotated.gtf")
#polyA_novel <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/02_novel_expressed/novel_expressed_polya.gtf")
#ribodepl_novel <- import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/02_novel_expressed/novel_expressed_ribodepleted.gtf")
#ref_gtf <-  import("/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/gene_annotations/UCSC/rheMac10_EBOV-Kikwit_UCSC.gtf")

# --------------------------------------------------------------------------------------
#                                   Merge Predictions 
# --------------------------------------------------------------------------------------

# -------------------------------------------
#                     CPC2 
# -------------------------------------------

# Load cpc2 prediction 
names(cpc2) <- c("ID","transcript_length ","peptide_length","Fickett_score" , " pI", "ORF_integrity", "coding_probability","label")
# Extract transcript ids 
cpc2$transcript_id <- unlist(lapply(as.character(cpc2$ID), function(x) strsplit(x,"[(]")[[1]][1]))
# Extract non coding predictions
pred_nc_ids <- cpc2[cpc2$label == "noncoding",]$transcript_id
cpc2_pre_lnc <- candidates[candidates$transcript_id %in% pred_nc_ids, ]

# -----------------------------------------------
# Remove uncorcondant pred
cpc2$gene_id <- unlist(lapply(as.character(cpc2$transcript_id), function(x) paste(unlist(strsplit(x, "[.]"))[1:2], collapse = ".")))

unconcordant_prediction <- cpc2 %>%  dplyr::group_by(gene_id) %>% dplyr::summarise(Unique_Elements =  dplyr::n_distinct(label)) %>%  dplyr::filter( Unique_Elements > 1)
cpc2_pre_lnc <- cpc2_pre_lnc[!(cpc2_pre_lnc$gene_id %in% unconcordant_prediction$gene_id)]
length(unique(cpc2_pre_lnc$gene_id))
# -----------------------------------------------


# -------------------------------------------
#                     CPAT 
# -------------------------------------------
# Compare which ones and how many overlap 
# Cutoff defined by the paper of CPAT as best for discerning coding and lnc
cutoff <-  0.364
cpat$label <- ifelse(cpat$Coding_prob > cutoff, "coding", "noncoding")
cpat$transcript_id <- unlist(lapply(as.character(cpat$ID), function(x) strsplit(x,"[(]")[[1]][1]))
# extract lncRNAs
cpat_nc <- cpat[cpat$label == "noncoding",]
cpat_pre_lnc <- candidates[candidates$transcript_id %in% cpat_nc$transcript_id, ]

# -----------------------------------------------
# Remove uncorcondant pred
cpat$gene_id <- unlist(lapply(as.character(cpat$ID), function(x) paste(unlist(strsplit(x, "[.]"))[1:2], collapse = ".")))
unconcordant_prediction <- cpat %>%  dplyr::group_by(gene_id) %>% dplyr::summarise(Unique_Elements =  dplyr::n_distinct(label)) %>%  dplyr::filter( Unique_Elements > 1)
cpat_pre_lnc <- cpat_pre_lnc[!(cpat_pre_lnc$gene_id %in% unconcordant_prediction$gene_id)]
length(unique(cpat_pre_lnc$gene_id))
# -----------------------------------------------



# -------------------------------------------
#                     CNIT 
# -------------------------------------------
# Compare which ones and how many overlap 
cnit$label <- cnit$index
cnit$transcript_id <- unlist(lapply(as.character(cnit$Transcript.ID), function(x) strsplit(x,"[(]")[[1]][1]))
cnit_nc <- cnit[cnit$label == "noncoding",]
cnit_pre_lnc <- candidates[candidates$transcript_id %in% cnit_nc$transcript_id, ]

# -----------------------------------------------
# Remove uncorcondant pred
cnit$gene_id <- unlist(lapply(as.character(cnit$Transcript.ID), function(x) paste(unlist(strsplit(x, "[.]"))[1:2], collapse = ".")))
cnit$transcript_id <- cnit$ID
unconcordant_prediction <- cnit %>%  dplyr::group_by(gene_id) %>% dplyr::summarise(Unique_Elements =  dplyr::n_distinct(label)) %>%  dplyr::filter( Unique_Elements > 1)
cnit_pre_lnc <- cnit_pre_lnc[!(cnit_pre_lnc$gene_id %in% unconcordant_prediction$gene_id)]
length(unique(cnit_pre_lnc$gene_id))
# -----------------------------------------------



intersection <- unique(intersect(intersect(cpat_pre_lnc$transcript_id, cpc2_pre_lnc$transcript_id), cnit_pre_lnc$transcript_id))
intersection_pre_lnc <- candidates[candidates$transcript_id %in% intersection, ]
ribodepl_novel <- intersection_pre_lnc

# --------------------------------------------------------------------------------------
#                                   Merge Predictions and reference
# --------------------------------------------------------------------------------------


# ---------------------------
#   Add Names 
# ---------------------------
add_gene_name <- function(gtf){
  missing_gene_name <- gtf[gtf$type == "exon"][is.na(gtf[gtf$type == "exon"]$gene_name)]                                
  transcripts <- gtf[gtf$type == "transcript"]
  # To the ones that have no gene name, find the gene name associated with the transcript ID and use it
  test <-  lapply(seq_along(missing_gene_name), function(x) transcripts[transcripts$transcript_id == missing_gene_name[x]$transcript_id]$gene_name)
  gtf[gtf$type == "exon"][is.na(gtf[gtf$type == "exon"]$gene_name)]$gene_name <- unlist(test)
  gtf[is.na(gtf$gene_name)]$gene_name <- paste(gtf[is.na(gtf$gene_name)]$gene_id, "-unknown", sep = "")
  return(gtf)
}
add_prefix <- function(gtf, prefix){
  gtf$gene_id <- paste(prefix, gtf$gene_id, sep = "") 
  gtf$transcript_id <- paste(prefix, gtf$transcript_id, sep = "") 
  return(gtf)
}

print("---0 ")

#gtf_ribodepl_newids <- add_prefix(ribodepl_novel, "ribodepl-")

# ---------------------------------------
# Merge Novel and Ref  together ( just append ) 
# ---------------------------------------
glst_ref_and_novel <-list(ref_gtf,ribodepl_novel)
gtf_ref_and_novel <- do.call(c, as(glst_ref_and_novel, "GRangesList"))


gtf_ref_and_novel_sort <- sort(gtf_ref_and_novel)
export(gtf_ref_and_novel_sort, outfile)  

print("---1 ")
mask <- is.na(gtf_ref_and_novel_sort$gene_name)
gtf_ref_and_novel_sort$gene_name
gtf_ref_and_novel_sort[mask]$gene_name <- paste(gtf_ref_and_novel_sort[mask]$gene_id, "-unknown", sep = "") 

print("---1 ")
if(is.null(gtf_ref_and_novel_sort$transcript_name)){
  gtf_ref_and_novel_sort$transcript_name <- NA
}
mask_t <- is.na(gtf_ref_and_novel_sort$transcript_name)
gtf_ref_and_novel_sort[mask_t]$transcript_name <-paste(gtf_ref_and_novel_sort[mask_t]$transcript_id, "-unknown", sep = "")  
export(gtf_ref_and_novel_sort, outfile_name)  


