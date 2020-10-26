#!/bin/sh

#input="/home/luisas/Desktop/cluster/gene_annotations/ensembl_release98/rheMac10/Macaca_mulatta.Mmul_10.98.gtf"
#input="/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all/03_novel_lncrnas/01_novel_expressed/novel_expressed_ribodepleted.gtf"
#input="/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all/03_novel_lncrnas/01_novel_expressed/novel_expressed_polyA.gtf"

name="rheMac10_EBOV_and_novel_genenames_exonid"

input="/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/$name.gtf"


human_reference="/home/luisas/Desktop/cluster/gene_annotations/ensembl_release100/homo_sapiens/Homo_sapiens.GRCh38.100.gtf"
#output_dir="/home/luisas/Desktop/cluster/data/01_Ebola-RNASeq_all/03_novel_lncrnas/02_final_catalogue/03_orthologs_novel/ribodepl"
output_dir="/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/04_orthologs/"

chain_rhemac_human="/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/gene_annotations/chains/rheMac10.hg38.rbest.chain"
chain_human_rhemac="/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/01_PreliminaryFiles_rheMac10/gene_annotations/chains/hg38.rheMac10.rbest.chain"

# Get bed with 4 fields for human reference also
output_dir_reference=$output_dir

mkdir -p $output_dir
# echo "convert human reference to bed6 "
Rscript ./Create_bed6.R $human_reference "$output_dir_reference/human_reference.bed6" "addbiotype"
bedtools sort -i "$output_dir_reference/human_reference.bed6" > "$output_dir_reference/human_reference.sorted.bed6"
rm "$output_dir_reference/human_reference.bed6"

echo "convert human reference to bed6 : DONE "

echo "convert input to bed6 "
# # Get bed with 4 fields for liftover
Rscript ./ExtractExons.R $input "$output_dir/$name.bed6" "lncrna"
bedtools sort -i "$output_dir/$name.bed6" > "$output_dir/$name.sorted.bed6"
rm "$output_dir/$name.bed6"
echo "convert input to bed6 : DONE "

#
# # Liftover macaque to human
liftOver -bedPlus=6 -tab "$output_dir/$name.sorted.bed6"  $chain_rhemac_human "$output_dir/$name.human.bed6" "unmapped"
bedtools sort -i "$output_dir/$name.human.bed6" > "$output_dir/$name.human.sorted.bed6"
rm "$output_dir/$name.human.bed6"
#
# # Liftover back to macaque
liftOver -bedPlus=6 -tab "$output_dir/$name.human.sorted.bed6"  $chain_human_rhemac "$output_dir/$name.human.macaque.bed6" "unmapped2"
bedtools sort -i "$output_dir/$name.human.macaque.bed6" > "$output_dir/$name.human.macaque.sorted.bed6"
rm "$output_dir/$name.human.macaque.bed6"
#
# # Intersect
bedtools intersect -wo -f 0.5 -r -s -a "$output_dir/$name.sorted.bed6" -b "$output_dir/$name.human.macaque.sorted.bed6" > "$output_dir/intersection.tab"

# Exploratory
# FindOverlaps of the human ones with human annotation so we can retrieve the gene biotype
# Left: human reference bed5 with gene biotype. Right: subset which was liftovered.
bedtools intersect -wo -f 0.5 -r -s -a  "$output_dir_reference/human_reference.sorted.bed6" -b "$output_dir/$name.human.sorted.bed6"> "$output_dir/intersection_human_macaquehuman.tab"



Rscript ./Create_bed6.R $human_reference "$output_dir_reference/human_reference_gene.bed6" "addbiotype"
bedtools sort -i "$output_dir_reference/human_reference_gene.bed6" > "$output_dir_reference/human_reference.gene.sorted.bed6"
rm "$output_dir_reference/human_reference_gene.bed6"

bedtools intersect -wo -f 0.5 -r -s -a  "$output_dir_reference/human_reference.gene.sorted.bed6" -b "$output_dir/$name.human.sorted.bed6"> "$output_dir/intersection_human_macaquehuman_genelevel.tab"
