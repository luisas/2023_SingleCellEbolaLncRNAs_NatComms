/*
*
* Nextflow pipeline for generation of preliminary files for the RNA-seq analys
* pipeline.
* Specifically it creates a merged assembly as well as annotation for the
* macaque and ebola virus.
* Moreover it generates the index necessary for running hisat2.
*
*/


// ------------ INPUT PARAMETERS -----------------------
params.prefix_rawdata = "/gpfs/projects/bsc83/Ebola/00_RawData/"
params.rhesus_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8/rheMac8.fa"
params.ebov_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.fa"
params.rheMac8_annotation_gz_file = "/gpfs/projects/bsc83/gene_annotation/ensembl_release97/rheMac8/Macaca_mulatta.Mmul_8.0.1.97.gtf.gz"
params.chromAlias = "/gpfs/projects/bsc83/Ebola/00_InformationFiles/gene_annotations/chromAlias.txt"

params.output_dir = "/gpfs/projects/bsc83/Ebola/data/"
params.prefix = "rheMac8_EBOV-Kikwit"
params.scripts="${baseDir}/scripts/"

//-------------- CREATE CHANNELS ------------------------
rhesus_genome_channel = Channel
                        .fromPath("${params.rhesus_genome}")
ebov_genome_channel = Channel
                      .fromPath("${params.ebov_genome}")
scripts=file("${params.scripts}")
rheMac8_annotation_gz_channel = Channel
                                  .fromPath("${params.rheMac8_annotation_gz_file}")
                                  .map{tuple(it.baseName.minus('.gtf'),it)}
chromAliasChannel = Channel
                    .fromPath("${params.chromAlias}")
ebov_annotation_channel = Channel
                          .fromPath("${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.gtf")

extract_ss_script = Channel
                    .fromPath("/gpfs/scratch/bsc83/bsc83768/bin/hisat2-2.1.0/hisat2_extract_splice_sites.py")
extract_exon_script = Channel
                      .fromPath("/gpfs/scratch/bsc83/bsc83768/bin/hisat2-2.1.0/hisat2_extract_exons.py")
gtfToGenePred_script_channel = Channel
                               .fromPath("/gpfs/scratch/bsc83/bsc83768/bin/ucsc_tools/gtfToGenePred")


/*
* Merge Assemblies of rheMac8 and ebov Virus to generate one merged assembly.
*/
process merge_assemblies {

    storeDir "${params.output_dir}/reference_assembly"

    input:
    file rheMac8 from rhesus_genome_channel
    file ebov from ebov_genome_channel

    output:
    file "${params.prefix}.fa" into (merged_assembly, merged_assembly_2,  merged_assembly_for_dictionary )

    script:
    """
    cat ${rheMac8} > ${params.prefix}.fa
    sed 's/KU182905.1/EBOV_Kikwit/g' ${ebov} >> ${params.prefix}.fa
    """

}

/*
* Create hisat2 indexes
*
* Build hisat2 indexes for joined rheMac8 and EBOV assemblies WITHOUT specifying
* exon and splice-sites files as done in the Tuxedo Protocol paper 2016
*
*/
process create_hisat2_indexes{

  storeDir "${params.output_dir}/indexes/hisat2"

  input:
  file assembly from merged_assembly

  output:
  file "${params.prefix}.*" into hisat2_indexes

  script:
  """
  hisat2-build-s ${assembly} ${params.prefix}
  """
}


/*
* Create dictionary for merged assembly with picard tools
*/
process create_dictionary{

  storeDir "${params.output_dir}/indexes/hisat2"

  input:
  file assembly from merged_assembly_for_dictionary

  output:
  file "*.dict" into dictionary_channel

  script:
  """
  java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar CreateSequenceDictionary R=${assembly}
  """
}


/*
* REQUIRED same chrom names in reference assembly and gene annotations
* Replace ensembl names with ucsc names in the rheMac8 annotation
*/
process modify_identifiers_annotation_rheMac{

  storeDir "${params.output_dir}/gene_annotations"

  input:
  set  annotation_name, file(rheMac8_annotation_gz) from rheMac8_annotation_gz_channel
  file ensembl_correspondence from ensembl_correspondence_channel
  file chromAlias from chromAliasChannel
  file scripts


  output:
  set file("rheMac8.ensembl_release97.gtf"), file("rheMac8.ensembl_release97.tab") into gene_annotations_rheMac_channel
  file "ensembl.ucsc.chrom_names_correspondance.tab" into ensembl_correspondence_channel


  script:
  """
  cat ${chromAlias}| awk '{ if (\$3 ~  /ensembl/) { print } }' | cut -f'1,2' > ensembl.ucsc.chrom_names_correspondance.tab
  # Replace ensembl names with ucsc names -> REQUIRED same chrom names in reference assembly and gene annotations
  awk 'FNR==NR {x2[\$1] = \$2; next} \$1 in x2 {OFS="\t";\$1=x2[\$1];print \$0}' ensembl.ucsc.chrom_names_correspondance.tab  <(zcat ${rheMac8_annotation_gz}) | cut -f1-8 > ${annotation_name}.column1-8.ucsc_chrom_names.txt
  # Create one file with description column
  zcat ${rheMac8_annotation_gz} | cut -f9 | awk '{if(\$1~/^#/){next}else{print}}' > ${annotation_name}.column9.txt
  # Complete parsing
  cat <(zcat  ${rheMac8_annotation_gz} | awk '{if(\$1~/^#/)print}') <( paste ${annotation_name}.column1-8.ucsc_chrom_names.txt ${annotation_name}.column9.txt) > rheMac8.ensembl_release97.gtf
  # Parse to tab delimited
  cat rheMac8.ensembl_release97.gtf |  ${scripts}/gtf2bed.py > rheMac8.ensembl_release97.tab
  """
}


/*
* Concatenate ensembl release 97 macaque annotation with EBOV gene annotation
*/
process merge_annotations{
  storeDir "${params.output_dir}/gene_annotations"

  input:
  file ebov_annotation from ebov_annotation_channel
  set file(rheMac_annotation), file(rheMac_bed) from gene_annotations_rheMac_channel

  output:
  file ("${params.prefix}.gtf") into (merged_annotation_channel, merged_annotation_channel_2)

  script:
  """
  cat  ${rheMac_annotation} ${ebov_annotation} > ${params.prefix}.gtf
  """

}

/*
* Extract exons and splice sites for annotated genes
*/
process extract_exons_ss{
 storeDir "${params.output_dir}/gene_annotations"

 input:
 file merged_annotation from merged_annotation_channel
 file extract_ss from extract_ss_script
 file extract_exon from extract_exon_script
 output:
 set file("${params.prefix}.exon.txt"), file("${params.prefix}.ss.txt") into extracted_exons_ss_channel

 shell:
 '''
 pyhton !{extract_ss} !{merged_annotation} > !{params.prefix}.ss.txt
 python !{extract_exon} !{merged_annotation} > !{params.prefix}.exon.txt
 '''
}



/*
* Create hisat2 indexes with exons and splicing sites information.
*/
process hisat2_ss_exon{

  storeDir "${params.output_dir}/indexes/hisat2_ss_exon"

  input:
  set file(extracted_ss), file(extracted_exons) from extracted_exons_ss_channel
  file merged_assembly from merged_assembly_2

  output:
  file "${params.prefix}.*" into hisat2_ss_exon_indexes

  script:
  """
  hisat2-build -p 8 --ss ${extracted_ss}  --exon ${extracted_exons} ${merged_assembly_2} ${params.prefix}
  """
}




// ------------------------------------
// ------------------------------------
// DO NOT DELETE YET!
// ------------------------------------
// ------------------------------------


// process parse_identifiers_correspondence{
//   storeDir "${params.output_dir}/gene_annotations"
//
//   input:
//   file chromAlias from chromAliasChannel
//
//   output:
//   file "ensembl.ucsc.chrom_names_correspondance.tab" into ensembl_correspondence_channel
//
//   script:
//   """
//   cat ${chromAlias}| awk '{ if (\$3 ~  /ensembl/) { print } }' | cut -f'1,2' > ensembl.ucsc.chrom_names_correspondance.tab
//   """
// }



// //# 4. Convert GFT to BED12
// //#/gpfs/scratch/bsc83/bsc83768/bin/ucsc_tools/gtfToGenePred rheMac8.ensembl_release97.gtf rheMac8.ensembl_release97.bed
// //#/gpfs/scratch/bsc83/bsc83768/bin/ucsc_tools/gtfToGenePred /gpfs/projects/bsc83/Ebola/data/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.gtf EBOV-Kikwit.bed
// ///gpfs/scratch/bsc83/bsc83768/bin/ucsc_tools/gtfToGenePred rheMac8_EBOV-Kikwit.gtf rheMac8_EBOV-Kikwit.tmp.bed
// //awk '{OFS="\t";print $2,$4,$5,$1,0,$3,$6,$7,0,$8,$9,$10}' rheMac8_EBOV-Kikwit.tmp.bed > rheMac8_EBOV-Kikwit.bed # This is the format required by infer_experiment.py from RNASeqQC
//
// process convert_gtf_to_bed12{
//   storeDir "${params.output_dir}/gene_annotations"
//
//   input:
//   file merged_assembly from merged_annotation_channel_2
//   file gtfToGenePred_script from gtfToGenePred_script_channel
//
//   output:
//   file "${params.prefix}.bed" into bed_channel
//
//   script:
//   """
//   bash ${gtfToGenePred_script} ${merged_assembly} ${params.prefix}.tmp.bed
//   awk '{OFS="\t";print \$2,\$4,\$5,\$1,0,\$3,\$6,\$7,0,\$8,\$9,\$10}' ${params.prefix}.tmp.bed > ${params.prefix}.bed
//   """
// }
