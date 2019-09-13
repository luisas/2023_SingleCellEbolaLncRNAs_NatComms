// Nextflow pipeline for data preparation  .
// possibly to be later merged with the analysis pipeline

params.prefix_rawdata = "/gpfs/projects/bsc83/Ebola/00_RawData/"
params.rhesus_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8/rheMac8.fa"
params.ebov_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.fa"

rhesus_genome_channel = Channel.fromPath("${params.rhesus_genome}")
ebov_genome_channel = Channel.fromPath("${params.ebov_genome}")


params.output_dir = "/gpfs/scratch/bsc83/bsc83024/test_pipeline"
/*
* Merge Assemblies (+ change name to match with annotation in ebov_genome)
* GTF gene annotation EBOV annotates as EBOV_Kikwit while in fasta file as
* Replace EBOV_Kikwit for chrom name in fasta file so it matches with annotation
*/
process merge_assemblies {

    storeDir "${params.output_dir}/reference_assembly"

    input:
    file rheMac8 from rhesus_genome_channel
    file ebov from ebov_genome_channel

    output:
    file "rheMac8_EBOV-Kikwit.fa" into (merged_assembly, merged_assembly_for_dictionary )

    script:
    """
    cat ${rheMac8} > rheMac8_EBOV-Kikwit.fa
    sed 's/KU182905.1/EBOV_Kikwit/g' ${ebov} >> rheMac8_EBOV-Kikwit.fa
    """

}



// /*
// * Create hisat2 indexes
// *
// * Build hisat2 indexes for joined rheMac8 and EBOV assemblies specifying
// * exon and splice-sites files as done in the Tuxedo Protocol paper 2016
// *
// */
// process create_hisat2_indexes{
//
//   storeDir "${params.output_dir}/indexes/hisat2"
//
//   input:
//   file assembly from merged_assembly
//
//   output:
//   file "*" into hisat2_indexes
//
//   script:
//   """
//   hisat2-build-s ${assembly} .
//   """
// }
//
//
// /*
// * Create dictionary for merged assembly
// */
// process create_dictionary{
//
//   storeDir "${params.output_dir}/indexes/hisat2"
//
//   input:
//   file assembly from merged_assembly_for_dictionary
//
//   output:
//   file "*.dict" into hisat2_indexes
//
//   script:
//   """
//   java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar CreateSequenceDictionary R=${assembly}
//   """
// }


/*
* REQUIRED same chrom names in reference assembly and gene annotations
* Replace ensembl names with ucsc names
*/
params.rheMac8_annotation_gz_file = "/gpfs/projects/bsc83/gene_annotation/ensembl_release97/rheMac8/Macaca_mulatta.Mmul_8.0.1.97.gtf.gz"
rheMac8_annotation_gz_channel = Channel
                                  .fromPath("${params.rheMac8_annotation_gz_file}")
                                  .map{it.baseName.minus(".gtf.gz"),it)))}
params.chromAlias = "/gpfs/projects/bsc83/Ebola/00_InformationFiles/gene_annotations/chromAlias.txt"
chromAliasChannel = Channel.fromPath("${params.chromAlias}")

process parse_identifiers_correspondence{
  storeDir "${params.output_dir}/gene_annotations"

  input:
  file chromAlias from chromAliasChannel

  output:
  file "ensembl.ucsc.chrom_names_correspondance.tab" into ensembl_correspondence_channel

  script:
  """
  #!/usr/bin/R
  chromAlias=read.table(${chromAlias})
  contigs <- chromAlias[chromAlias$V3=="refseq",2]
  x1=chromAlias[sapply(chromAlias$V3, function(i) "ensembl" %in% i),]
  x2=chromAlias[sapply(chromAlias$V3, function(i) "genbank,ensembl" %in% i),]
  z=rbind.data.frame(x1,x2)[,c(1:2)]
  write.table(z,"ensembl.ucsc.chrom_names_correspondance.tab",col.names=F, row.names=F,sep="\t",quote=F)
  """
}


process gene_annotation{

  storeDir "${params.output_dir}/gene_annotations"

  input:
  set  annotation_name, file(rheMac8_annotation_gz) from rheMac8_annotation_gz_channel
  file ensembl_correspondence from ensembl_correspondence_channel

  output:
  set file("rheMac8.ensembl_release97.gtf"), file("rheMac8.ensembl_release97.tab") into gene_annotations

  script:
  """
  # Replace ensembl names with ucsc names -> REQUIRED same chrom names in reference assembly and gene annotations
  awk 'FNR==NR {x2[$1] = $2; next} $1 in x2 {OFS="\t";$1=x2[$1];print $0}' ${ensembl_correspondence}  <(zcat ${rheMac8_annotation_gz}) | cut -f1-8 > ${annotation_name}.column1-8.ucsc_chrom_names.txt
  # Create one file with description column
  zcat ${rheMac8_annotation_gz} | cut -f9 | awk '{if($1~/^#/){next}else{print}}' > ${annotation_name}.column9.txt
  # Complete parsing
  cat <(zcat  ${rheMac8_annotation_gz} | awk '{if($1~/^#/)print}') <( paste ${annotation_name}.column1-8.ucsc_chrom_names.txt ${annotation_name}.column9.txt) > rheMac8.ensembl_release97.gtf
  # Parse to tab delimited
  cat rheMac8.ensembl_release97.gtf | ./gtf2bed.py > rheMac8.ensembl_release97.tab
  """
}


// missing
// # 2. Concatenate ensembl release 97 macaque annotation with EBOV gene annotation
// cat  rheMac8.ensembl_release97.gtf /gpfs/projects/bsc83/Ebola/data/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.gtf > rheMac8_EBOV-Kikwit.gtf

process merge_annotations{


}

// # 3. Extract exons and splice sites for annotated genes
// sbatch extract_exons_ss.sh
process extract_exons_ss{

}

//# 4. Convert GFT to BED12
//#/gpfs/scratch/bsc83/bsc83768/bin/ucsc_tools/gtfToGenePred rheMac8.ensembl_release97.gtf rheMac8.ensembl_release97.bed
//#/gpfs/scratch/bsc83/bsc83768/bin/ucsc_tools/gtfToGenePred /gpfs/projects/bsc83/Ebola/data/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.gtf EBOV-Kikwit.bed
///gpfs/scratch/bsc83/bsc83768/bin/ucsc_tools/gtfToGenePred rheMac8_EBOV-Kikwit.gtf rheMac8_EBOV-Kikwit.tmp.bed
//awk '{OFS="\t";print $2,$4,$5,$1,0,$3,$6,$7,0,$8,$9,$10}' rheMac8_EBOV-Kikwit.tmp.bed > rheMac8_EBOV-Kikwit.bed # This is the format required by infer_experiment.py from RNASeqQC

process convert_gtf_to_bed{

}


// /*
// * to be done when annotations are ready
// */
// process hisat2_ss_exon{
//
//   storeDir "${params.output_dir}/indexes/hisat2"
//
//   input:
//   file assembly from merged_assembly_for_dictionary
//
//   output:
//   file "*.dict" into hisat2_indexes
//
//   script:
//   """
//   hisat2-build -p 8 --ss /gpfs/projects/bsc83/Ebola/indexes/hisat2_ss_exon/rheMac8_EBOV-Kikwit.ss.txt  --exon /gpfs/projects/bsc83/Ebola/indexes/hisat2_ss_exon/rheMac8_EBOV-Kikwit.exon.txt /gpfs/projects/bsc83/Ebola/indexes/hisat2_ss_exon/rheMac8_EBOV-Kikwit.fa /gpfs/projects/bsc83/Ebola/indexes/hisat2_ss_exon/rheMac8_EBOV-Kikwi
//   """
// }
