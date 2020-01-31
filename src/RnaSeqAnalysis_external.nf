#!/usr/bin/env nextflow

/*
 * TODO COPYRIGHT
 */

 log.info "=================================================================================="
 log.info "     RNASeq analysis pipeline & de Novo Transcriptome assembly"
 log.info "=================================================================================="

 /* --------------------------------------
  *  INPUT PARAMETERS defaults
  * --------------------------------------
  */

// ------------- RAW DATA

// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirDataAnnot = "/gpfs/projects/bsc83/Data"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"

params.dataset_fastq_dir_ext = "${params.dirData}/00_RawData/extrenal_rhesus_RNA-Seq/SRP016501-rhesus/01_fastq/*/*/*/"
params.dataset_fastq_dir_ext_two = "${params.dirData}/00_RawData/extrenal_rhesus_RNA-Seq/nhprt-wholeblood-chineserhesus/01_fastq/*/*/*/"

// Folder where the output directories of the pipeline will be placed
params.output_dir = "${params.dirData}/01_Ebola-RNASeq/02_RNA-Seq_external_correctstrand/"
//params.output_dir = "/gpfs/projects/bsc83/Projects/Ebola/test_data/04_test/"

// ---------- PRELIMINARY FILES

params.rehmacrelease = "rheMac10"
params.preliminary_files_dir="${params.dirData}/01_Ebola-RNASeq/01_PreliminaryFiles_${params.rehmacrelease}"
params.assembly_name = "${params.rehmacrelease}_EBOV-Kikwit"

// This ones should be automatically retrieved but can be changed in need
// Assemblies and indexes
params.reference_assembly = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.fa"
params.reference_assembly_fai = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.fa.fai"
params.dict = "${params.preliminary_files_dir}/reference_assembly/${params.assembly_name}.dict"
params.hisat2_indexes = "${params.preliminary_files_dir}/indexes/hisat2/${params.assembly_name}"
params.star_indexes = "${params.preliminary_files_dir}/indexes/star/"

params.ss = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.ss.txt"
// Gene annotation
params.bed = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.bed"
params.gtf = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.gtf"
params.gtf_ref_merged =  "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}_and_novel.gtf"
//params.bed_rheMac = "${params.dirData}/00_RawData/pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/rheMac8/Ensembl/rheMac8.Ensembl.bed"
params.bed_rheMac = "${params.preliminary_files_dir}/gene_annotations/${params.assembly_name}.bed12"


// If set to true the quality assessment will be computed with fastqc
params.fastqc=true


// -----------------------------------------------

/* Channel for all the unmapped bam files present in any of the dataset_bam_dir
* subdirectories.
* In order, it maps:
* - lane
* - lane_number
* - dataset_name
* - dayPostInfection
* - tissue
* - sample
* - complete_id (${dataset_name}_${dayPostInfection}_${tissue}_${sample}_l${lane_number})
* - the file itself.
*/


// Fastq1 have no lanes
fastq_one = Channel.fromFilePairs("${params.dataset_fastq_dir_ext}/*_{1,2}.fastq.gz")
              .ifEmpty('fastq files directory is empty')
              .map { tuple(it[1][0].toString().split('/').reverse()[5],
                           it[1][0].toString().split('/').reverse()[3],
                           it[1][0].toString().split('/').reverse()[2],
                           it[1][0].toString().split('/').reverse()[1],
                           "1",
                           it[1][0].toString().split('/').reverse()[5] + "_" +it[1][0].toString().split('/').reverse()[3]+ "_" +it[1][0].toString().split('/').reverse()[2]+ "_" +it[1][0].toString().split('/').reverse()[1]+ "_l1",
                           it[1] ) }


fastq_two = Channel.fromFilePairs("${params.dataset_fastq_dir_ext_two}/*_R{1,2}_001.fastq.gz")
              .ifEmpty('fastq files directory is empty')
              .map { tuple(it[1][0].toString().split('/').reverse()[5],
                           it[1][0].toString().split('/').reverse()[3],
                           it[1][0].toString().split('/').reverse()[2],
                           it[1][0].toString().split('/').reverse()[1],
                           it[1][0].baseName.split('_')[4].split('00')[1],
                           it[1][0].toString().split('/').reverse()[5] + "_" +it[1][0].toString().split('/').reverse()[3]+ "_" +it[1][0].toString().split('/').reverse()[2]+ "_" +it[1][0].toString().split('/').reverse()[1]+ "_l"+it[1][0].baseName.split('_')[4].split('00')[1],
                           it[1] ) }


fastqs = fastq_one.mix(fastq_two)
fastqs.into{ dataset_fastq_ext; fastq_files_for_mapping; printing }



// ------------ CHANNELS Creation
Channel.fromPath("${params.hisat2_indexes}*").into{ indexesForMapping; indexesForMapping2 }
Channel.fromPath("${params.star_indexes}").set{ star_indexes}
Channel.fromPath("${params.reference_assembly}").into{ reference_assembly_channel; reference_genome_ch;reference_genome_ch_2; reference_genome_ch_3}
Channel.fromPath("${params.reference_assembly_fai}").set{reference_genome_fai_ch}
Channel.fromPath("${params.ss}").set{ known_ss; }
Channel.fromPath("${params.bed}")
                                 .into{ annotation_bed; annotation_bed2 }
Channel.fromPath("${params.bed_rheMac}")
                                 .set{ annotation_bed_rhemac }
Channel.fromPath("${params.dict}").set{dictionary_channel}
Channel.fromPath("${params.gtf}").into{ gtfChannel1; gtfChannel2; gtfChannel3; gtfChannel4; gtfChannel5; gtfChannel6; gtfChannel7}
Channel.fromPath("${params.gtf_ref_merged}").set{ gtfANDnovel }

/*  ----------------------------------------------------------------------
*   ----------------------------------------------------------------------
*                       BEGINNING OF THE PIPELINE
*   ----------------------------------------------------------------------
*   ----------------------------------------------------------------------
*/



/*  java -Xmx4g  -Djava.io.tmpdir=$TMPDIR -jar /apps/PICARD/2.20.0/picard.jar \
* STEP 2: Generate quality assessment with fastqc.
*/
process generate_fastqc{

  tag "${complete_id}"
  storeDir "${params.output_dir}/02_fastqc/$dataset_name/$tissue/$dayPostInfection/$sample"


  input:
  set  dataset_name, tissue,dayPostInfection, sample, lane, complete_id,
      file(fastq_pair) from dataset_fastq_ext

  output:
  file "*" into fastqcs

  script:
  """
  fastqc ${fastq_pair}  --extract
  """

}

/*
* STEP 3.1: Mapping with HISAT.
* It maps the reads of the generated fastq files
* to the assembly defined in params.assembly
* The indexes were computed WITHOUT the help of known splice sites
* the known splice sites are now used as input for the mapper.
*/
process mapping_hisat{

  label 'big_mem'
  cpus 8
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  file indexes from indexesForMapping.collect()
  file ss from known_ss.collect()
  set dataset_name, tissue, dayPostInfection, sample, lane, complete_id,
      file(fastq_pair) from fastq_files_for_mapping

  output:
  set dataset_name, tissue,dayPostInfection, sample, lane, complete_id,
      file("${complete_id}.novel_ss.txt"), file("${complete_id}.hisat2_summary.txt"),
      file("${complete_id}.sam") into mapped_sam

  /*
  *  --known-splicesite-infile <path>   provide a list of known splice sites
  *  --novel-splicesite-outfile <path>  report a list of splice sites
  *  --dta                              reports alignments tailored for transcript assemblers
  *  --summary-file                     Print alignment summary to this file.
  */

  script:
  """
  hisat2 -p ${task.cpus} -x ${params.assembly_name} -1  ${fastq_pair[0]} \
                      -2 ${fastq_pair[1]} --known-splicesite-infile ${ss} \
                      --novel-splicesite-outfile ${complete_id}.novel_ss.txt \
                      --downstream-transcriptome-assembly \
                      --time --summary-file ${complete_id}.hisat2_summary.txt \
                      --rna-strandness RF > ${complete_id}.sam
  """

}


/*
* STEP 3.2: Sort and index the mapped bam files
*
* The sorting is done with picard as we encountered some problems when running
* samtools and afterwards other picardTools on top of the sorted bam.
*/
process sort_bam{

  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"


  input:
  set dataset_name, tissue,dayPostInfection, sample, lane, complete_id,
      file(novel_ss),file(summary),
      file(sam) from mapped_sam

  output:
  set file("${complete_id}.bam")   into sorted_and_index_bam

  script:
  """
  samtools sort ${sam} -O bam -o ${complete_id}.bam
  """
}

// // In this channel we group the files by the complete_id after removing the lane.
// // We have now a set of bams per sample from all the lanes.
sorted_and_index_bam.map { file ->
        def key = file.name.toString().tokenize('_').init()
        key_string = key.join("_")
        key << key_string
        return tuple(key, file)
     }
     .groupTuple()
     .map{ tuple(it.get(0).get(0),
                 it.get(0).get(1),
                 it.get(0).get(2),
                 it.get(0).get(3),
                 it.get(0).get(4),
                 it.get(1))}
    .set{ groups_lanes }


process merge_lanes{

    tag "${complete_id}"
    storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

    input:
    set dataset_name, tissue, dayPostInfection, sample, complete_id, file(samples) from groups_lanes

    output:
    set dataset_name, dayPostInfection, tissue, sample, complete_id, file("${complete_id}.UMI.bam"), file("${complete_id}.UMI.bam.bai") into (merged_bylanes, raw_bams)

    script:
    """
    samtools merge -f ${complete_id}.UMI.bam  ${samples}
    # Merge mapped bam files with UMI info
    samtools index  ${complete_id}.UMI.bam
    """
}


/*
*   STEP 6 filter bams
*   Filter low-quality reads and retain properly-paired uniquely mapped reads
*   f3 : read paired & read mapped in proper pair
*   q60: 60 should be uniquely mapped reads have MAPQ of 60.
*   Now it is a stringent filtering. We may wanna make it looser for the lncRNA
*   annotation pipeline.
*/
process filter_bams_samtools{

  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(bam), file(bai) from merged_bylanes

  output:
  set dataset_name, tissue,dayPostInfection, sample, complete_id, file("${complete_id}.f3.q60.bam"), file("${complete_id}.f3.q60.bam.bai") into (filtered_merged_bams_1, filtered_merged_bams_2,filtered_merged_bams_3, filtered_merged_bams_4, filtered_merged_bams_5, filtered_merged_bams_6)

  script:
  """
  samtools view -b -f3 -q 60  ${bam} > ${complete_id}.f3.q60.bam
  samtools index ${complete_id}.f3.q60.bam
  """
}


process runRSeQC{

  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/rseqc"

  input:
  set dataset_name, tissue,dayPostInfection, sample, complete_id,
      file(bam),
      file(bai) from filtered_merged_bams_6
  file(bed) from  annotation_bed_rhemac.collect()

  output:
  file "*" into read_distribution_channel

  script:
  """
  read_distribution.py -i ${bam} -r ${bed} > ${complete_id}.UMI.f3.q60.read_distribution.txt
  infer_experiment.py -i ${bam} -r ${bed} > ${complete_id}.infer_experiment.txt
  junction_annotation.py -i ${bam} -o ${complete_id}.rseqc -r ${bed}
  bam_stat.py -i ${bam} > ${complete_id}.bam_stat.txt
  junction_saturation.py -i ${bam} -o ${complete_id}.rseqc -r ${bed} > ${complete_id}.junction_annotation_log.txt
  inner_distance.py -i ${bam} -o ${complete_id}.rseqc -r ${bed}
  read_duplication.py -i ${bam} -o ${complete_id}.read_duplication
  """
}


/*
*  DE NOVO TRANSCRIPTOME ASSEMBLY
*/


// -A <gene_abund.tab>	Gene abundances will be reported (tab delimited format) in the output file with the given name.
// -B	This switch enables the output of Ballgown input table files (*.ctab) containing coverage data for the reference transcripts given with the -G option.
// --fr	Assumes a stranded library fr-secondstrand.
// -f <0.0-1.0>	Sets the minimum isoform abundance of the predicted transcripts as a fraction of the most abundant transcript assembled at a given locus.

params.transcriptome_assembly = true
process DeNovoAssembly{

  cpus 4
  tag "${complete_id}"
  storeDir "${params.output_dir}/04_stringtie/$dataset_name/$tissue/$dayPostInfection/$sample"

  when:
  params.transcriptome_assembly

  input:
  set dataset_name, tissue,dayPostInfection, sample, complete_id, file(bam), file(bai) from filtered_merged_bams_1
  file gtf from gtfChannel4.collect()

  output:
  val dataset_name into dataset_name_ch
  file "*.stringtie.gtf" into stringTie_channel

  script:
  file_prefix = get_file_name_no_extension(bam.name)
  """
  stringtie ${bam} -G ${gtf} -o ${file_prefix}.stringtie.gtf --rf -p ${task.cpus}
  """
}

/*
* Merge all the assemblies Reference Guided
*/
process StringTie_Merge_Reference_Guided{

  cpus 48
  storeDir "${params.output_dir}/04_stringtie/"

  when:
  params.transcriptome_assembly

  input:
//  val dataset_name from dataset_name_ch
  file(stringtie_gtfs) from stringTie_channel.collect()
  file reference_gtf from gtfChannel5

  output:
  file "stringtie_merged_reference_guided.gtf" into (merged_denovo_assmebly, merged_de_novo_assembly_2)
  //val dataset_name into dataset_name_ch2

  script:
  """
  stringtie --merge -p ${task.cpus} -F 1.0 -o stringtie_merged_reference_guided.gtf -G ${reference_gtf} ${stringtie_gtfs}
  """
}


// /*
// * Obtain stats over StringTie output
// */
// process gffCompare2{
//
//   storeDir "${params.output_dir}/04b_gffCompare"
//
//   input:
//   file merged_gtf from merged_denovo_assmebly
//   file reference_gtf from gtfChannel6
// //  val dataset_name from dataset_name_ch2
//
//   output:
//   file("merged*") into gff_compare_output_channel2
//
//   script:
//   """
//   gffcompare -R -r ${reference_gtf} -o merged ${merged_gtf}
//   """
// }


merged_de_novo_assembly_2.into{ merged_de_novo_assembly_3; merged_de_novo_assembly_4; merged_de_novo_assembly_5}


//params.data_folder = "/Users/luisasantus/Desktop/mn_cluster/mount_dirs/projects/data"
params.proj_data_folder = "/gpfs/projects/bsc83/Data"
params.dir_ebola_data = "${params.proj_data_folder}/Ebola"
params.data_folder = "${params.dir_ebola_data}/01_Ebola-RNASeq"
//params.assembly_name = "rheMac10_EBOV-Kikwit"
//params.output_dir = "${params.dir_ebola_data}/01_Ebola-RNASeq/02_RNA-Seq_rheMac10/"


params.ref_gtf = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/${params.assembly_name}.gtf"
//params.merged_gtf = "${params.data_folder}/02_RNA-Seq_rheMac10/04_stringtie/stringtie_merged_reference_guided.gtf"

// TRAINING DATA
// Trained with human data!
params.known_mrna = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownMrnas_nochr.gtf"
params.known_lncrna = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownlncRNAs_nochr.gtf"

params.rhesus_genome = "${params.dir_ebola_data}/01_Ebola-RNASeq/01_PreliminaryFiles_rheMac10/reference_assembly/rheMac10_EBOV-Kikwit.fa"

Channel.fromPath("${params.rhesus_genome}").into{ fasta_reference_channel; fasta_reference_channel2}
Channel.fromPath("${params.merged_gtf}").into{ merged_gtf_channel; merged_gtf_channel_2 }
Channel.fromPath("${params.ref_gtf}").into{ ref_gtf_channel; ref_gtf_channel_2 }
Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }


log.info "=============================================="
log.info "            lncRNA annotation"
log.info "=============================================="

process feelnc_filter{

  storeDir "${params.output_dir}/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc"

  input:
  file merged_gtf from merged_de_novo_assembly_3
  file reference_gtf from ref_gtf_channel

  output:
  file("candidate_lncRNA.gtf") into candidates

  script:
  """
  FEELnc_filter.pl -i ${merged_gtf} \
                   -a ${reference_gtf} \
                   -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
  """
}

process feelnc_codpot{
  cpus 48
  storeDir "${params.output_dir}/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc"

  input:
  file candidate_lncrna from candidates
  file known_mrna from known_mrna_channel
  file known_lncrna from known_lncrna_channel
  file reference_genome from fasta_reference_channel

  output:
  set file("rheMac10_EBOV-Kikwit.fa.index"), file("feelnc_codpot_out") into coding_potentials

  script:
  """
  FEELnc_codpot.pl -i ${candidate_lncrna} -a ${known_mrna} -l ${known_lncrna} \
                   -g ${reference_genome} --proc ${task.cpus}
  """

}


process gffCompare{

  storeDir "${params.output_dir}/05_lncrnaAnnotation_no_l_option/feelnc_gencode_linc/01_gffcompare"

  input:
  set index, codpot from coding_potentials
  file reference_gtf from ref_gtf_channel_2

  output:
  file("merged*") into gff_compare_output_channel

  script:
  """
  gffcompare -R -r ${reference_gtf} -o merged ${codpot}/candidate_lncRNA.gtf.lncRNA.gtf
  """
}


// process getHTseqCountMerged{
//
//   storeDir "${params.output_dir}/06_counts/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"
//
//   input:
//   file gtf from merged_de_novo_assembly_4.collect()
//   set dataset_name, dayPostInfection, tissue, sample, lane, complete_id,
//       file(bam),
//       file(bai) from filtered_merged_bams_2
//
//   output:
//   file "${bam_prefix}.HTseq.gene_counts.tab" into htseqCountsMerged_channel
//
//   script:
//   bam_prefix = get_file_name_no_extension(bam.name)
//   """
//   htseq-count -f bam -r pos -s yes -t exon -i gene_id ${bam} ${gtf} > ${bam_prefix}.HTseq.gene_counts.tab
//   """
// }




// process getHTseqCountRef{
//
//   storeDir "${params.output_dir}/06_counts_ref/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"
//
//   input:
//   file gtf from gtfChannel1.collect()
//   set dataset_name, dayPostInfection, tissue, sample, lane, complete_id,
//       file(bam),
//       file(bai) from filtered_merged_bams_3
//
//   output:
//   file "${bam_prefix}.HTseq.gene_counts.tab" into htseqCountsRef_channel
//
//   script:
//   bam_prefix = get_file_name_no_extension(bam.name)
//   """
//   # Calc the counts for the umi_dedup
//   htseq-count -f bam -r pos -s yes -t exon -i gene_id ${bam} ${gtf} > ${bam_prefix}.HTseq.gene_counts.tab
//   """
// }


// I merge the novel with the reference

// process getHTseqCountALL{
//
//   cpus 48
//   storeDir "${params.output_dir}/06_counts_all/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"
//
//   input:
//   file gtf from gtfANDnovel.collect()
//   set dataset_name, dayPostInfection, tissue, sample, lane, complete_id,
//       file(bam),
//       file(bai) from filtered_merged_bams_4
//
//   output:
//   file "${bam_prefix}.HTseq.gene_counts.tab" into htseqCountsALL_channel
//
//   script:
//   bam_prefix = get_file_name_no_extension(bam.name)
//   """
//   htseq-count -f bam -s yes -t exon -i gene_id ${bam} ${gtf} > ${bam_prefix}.HTseq.gene_counts.tab
//   """
// }
// //






//params.fix_merged = "${params.output_dir}/04_stringtie/Zyagen/Zyagen_stringtie_merged_reference_guided.gtf"
//Channel.fromPath("${params.fix_merged}").set{ fixed_zyagen_merged_stringtie}

/*
* QUANTIFICATION STEP - Stringtie abundances
*/
// process StringTie_abundances_umis{
//
//   cpus 8
//   storeDir "${params.output_dir}/06_abundances/$dataset_name/$tissue/$dayPostInfection/$sample/umis"
//
//   input:
//   //file stringtie_merged from fixed_zyagen_merged_stringtie.collect()
//   file stringtie_merged from merged_de_novo_assembly_5.collect()
//   set dataset_name, dayPostInfection, tissue, sample, lane, complete_id,
//       file(bam),
//       file(bai) from dedup_umis_3
//
//   output:
//   set complete_id, file("${file_prefix}_abundance.gtf") into abundances_stringtie_umis
//   file "*.ctab" into ctabs_umis
//   val dataset_name into dataset_name_channel
//
//   script:
//   file_prefix = get_file_name_no_extension(bam.name)
//   """
//   stringtie -e -B -p ${task.cpus} -G ${stringtie_merged} -A ${file_prefix}.ctab -o ${file_prefix}_abundance.gtf ${bam}
//   """
// }
//
//
//
// /*
// * CONVERT : get gene count matrix from StringTie output
// */
// params.script_prepde="${baseDir}/scripts/prepDE.py"
// Channel.fromPath("${params.script_prepde}").set{ prepdescript }
// abundances_stringtie_umis.into{abundances_stringtie_umis1 ; abundances_stringtie_umis2}
//
//
// process prepareCountMatrices{
//
//   storeDir "${params.output_dir}/06_quantification_stringtie"
//
//   input:
//   file fileList from abundances_stringtie_umis1.collect()
//   file samples_list from abundances_stringtie_umis2.flatten().collate(2).collectFile(name: 'sample.txt', newLine: true){it[0] +"\t"+ it[1].baseName+".gtf"}
//   file script from prepdescript.collect()
//
//   output:
//   //file "*" into out_prepde_umis
//   file "*"
//
//   script:
//   """
//   python2 ${script} -i ${samples_list}
//   """
// }


workflow.onComplete {
	println ( workflow.success ? "Done!" : "Oops .. something went wrong" )
}

/*   -------------------------------
*           Groovy Functions
*    -------------------------------
*/

def remove_lane_from_id(String id){
  return id.split("_").init().join("_")
}

def get_file_name_no_extension(String filename){
  return filename.split("\\.").init().join('.')
}
