library(rtracklayer)
library(GenomicRanges)

reference_ucsc <- import("/home/luisas/Desktop/cluster/data/99_BroadAnnotation/01_stringtie_assembly_merged/01_gffCompare/merged.annotated.gtf")


transcripts <-reference_ucsc[reference_ucsc$type == "transcript", ]

length(unique(transcripts[transcripts$class_code %in% c("u", "x"), ]$transcript_id))


candidates <- import("/home/luisas/Desktop/cluster/data/99_BroadAnnotation/02_FEELnc_prediction_old//candidate_lncRNA.gtf")

length(unique(candidates$transcript_id))
