

gtf_macaque="/home/luisas/Desktop/cluster/gene_annotations/ensembl_release98/rheMac10/Macaca_mulatta.Mmul_10.98.gtf"

output_dir="/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all/99_MitochondrialInvestigation"

chain="/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all/01_PreliminaryFiles_rheMac10/gene_annotations/chains/rheMac10ToRheMac8.over.chain"

convert_script="/home/luisas/Desktop/cluster/proj/code/ebola/src/scripts/convert_annotation.R"
gtf_human="$output_dir/Macaca_mulatta.Mmul_8.0.1.92.gtf"
# # Convert human seqLevel to UCSC so we can use the chains
Rscript "$convert_script" "$gtf_macaque" "UCSC" "$output_dir/macaque_ucsc.gtf"
#
# # Liftover Macaque to human
liftOver -gff "$output_dir/macaque_ucsc.gtf"  $chain "$output_dir/rhemac10to8.gff" "unmapped"

# Extract human mitochondrial annotation
Rscript ./extractMitochondrial.R "$gtf_human" "$output_dir/rhemac8_mt.gtf"

# Intersect human mitochondrial annotation and macaque
bedtools intersect -wo -a "$output_dir/rhemac8_mt.gtf" -b "$output_dir/rhemac10to8.gff" > "$output_dir/intersection_new.tab"
