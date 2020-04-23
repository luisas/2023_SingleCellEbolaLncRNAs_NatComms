/*
*
* Nextflow pipeline for generation of preliminary files for the RNA-seq analys
* pipeline.
* Specifically it creates a merged assembly as well as annotation for the
* macaque and ebola virus.
* Moreover it generates the index necessary for running hisat2.
*
*/
log.info "=============================================="
log.info " RNASeq analysis and annotation pipeline "
log.info "=============================================="

// ------------ INPUT PARAMETERS -----------------------

params.macaque_assembly_name = "rheMac10"
params.prefix = "${params.macaque_assembly_name}_EBOV-Kikwit"
params.release = "release98"

// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
//params.dirData = "/home/luisas/Desktop/cluster/data"
params.dataset_bam_dir_zyagen = "${params.dirData}/00_RawData/pardis_shared_data/sabeti-txnomics/alin/190713_Zyagen-longRNA/tmp/00_demux/bams_per_lane/*/"
params.dataset_bam_dir_batch = "${params.dirData}/00_RawData_links/"

params.prefix_data = "/gpfs/projects/bsc83/Data/"
//params.prefix_data = "/home/luisas/Desktop/cluster/"
params.rhesus_gtf = "${params.prefix_data}gene_annotation/ensembl_${params.release}/${params.macaque_assembly_name}/*.gtf"
params.rhesus_genome = "${params.prefix_data}assemblies/ensembl/release-98/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa"



params.prefix_rawdata = "${params.dirData}/00_RawData/"
params.ebov_genome = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.fa"
params.ebov_gtf = "${params.prefix_rawdata}pardis_shared_data/sabeti-txnomics/shared-resources/HISAT2/EBOV-Kikwit/KU182905.1.gtf"


// TRAINING DATA
params.data_folder = "${params.dirData}/01_Ebola-RNASeq"
params.known_mrna = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownMrnas_nochr.gtf"
params.known_lncrna = "${params.data_folder}/01_PreliminaryFiles_rheMac10/gene_annotations/gencode.v26.GRCh38.annotation_knownlncRNAs_nochr.gtf"
Channel.fromPath("${params.known_mrna}").set{ known_mrna_channel}
Channel.fromPath("${params.known_lncrna}").set{ known_lncrna_channel }


// Here it should not be really changed
params.output_dir_preliminary = "${params.dirData}/01_Ebola-RNASeq_all/01_PreliminaryFiles_${params.macaque_assembly_name}/"
params.output_dir_name = "02_RNA-Seq_external"
params.output_dir = "${params.dirData}/01_Ebola-RNASeq_all/${params.output_dir_name}/"


params.scripts="${baseDir}/scripts/"

params.umi = "true"
params.strandness = "FR"


if( "${params.strandness}" == "FR" ){
  params.htseqsense="yes"
}
if( "${params.strandness}" == "RF" ){
  params.htseqsense="reverse"
}



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

//  STILL HARDCODED :/ !!
// I would need to rename all the bams of Batch and Zyagen
// either rename the files or hardcode channels
// ---------------------------------

if( "${params.umi}" == "true" ){
  unmapped_bams_1 = Channel
                  .fromPath("${params.dataset_bam_dir_zyagen}/*_long.bam")
                  .ifEmpty('bam files directory is empty')
                  .map{ tuple(it.baseName.split('_')[0] + "_" + it.baseName.split('_')[2].split('-')[0] +"_"+ it.baseName.split('_')[1] +"_"+ it.baseName.split('_')[2].split('-')[1] +"_"+ "l" +  it.parent.name.split('\\.')[1],
                              it)}

  unmapped_bams_2 = Channel
                  .fromPath("${params.dataset_bam_dir_batch}/*.bam")
                  .ifEmpty('bam files directory is empty')
                  .map{ tuple(it.baseName.split('_')[0] + "_" + it.baseName.split('_')[1]+ "_" +it.baseName.split('_')[2] + "_" +it.baseName.split('_')[3] + "_" + "l"+it.baseName.split('_')[4].split('\\.')[1],
                              it)}

  unmapped_bams = unmapped_bams_1.mix(unmapped_bams_2)
}
// ---------------------------------

params.fastqs  = "/gpfs/projects/bsc83/Data/Ebola/00_RawData/extrenal_rhesus_RNA-Seq/*/*/Brain/*/*/*.{1,2}.fastq.gz"
fastqs = Channel.fromFilePairs("${params.fastqs}")
                .ifEmpty("No fastqs found")
                .map { tuple(it[0].split('_')[0],
                             it[0].split('_')[1],
                             it[0].split('_')[2],
                             it[0].split('_')[3],
                             it[0].split('_')[4].split("l")[1],
                             it[0].split('_')[0] + "_" +\
                             it[0].split('_')[1] + "_" +\
                             it[0].split('_')[2] + "_" +\
                             it[0].split('_')[3] + "_" +\
                             it[0].split('_')[4],
                             it[1] ) }


fastqs.into{ fastq_files_for_qc; fastq_files_for_mapping; printing }
printing.subscribe{ println "$it" }


//-------------- CREATE CHANNELS ------------------------
rhesus_genome_channel = Channel
                        .fromPath("${params.rhesus_genome}")
ebov_genome_channel = Channel
                      .fromPath("${params.ebov_genome}")
scripts=file("${params.scripts}")
rheMac_annotation_channel = Channel
                                  .fromPath("${params.rhesus_gtf}")

ebov_annotation_channel = Channel
                          .fromPath("${params.ebov_gtf}")

extract_ss_script = Channel
                    .fromPath("${params.scripts}/hisat2_extract_splice_sites.py")
extract_exon_script = Channel
                    .fromPath("${params.scripts}/hisat2_extract_exons.py")

gtfToGenePred_script_ch = Channel
                      .fromPath("${params.scripts}/gtfToGenePred")
genePredToBed_script_ch = Channel
                      .fromPath("${params.scripts}/genePredToBed")




log.info "=============================================="
log.info " Preliminary Files preparation  "
log.info "=============================================="
/*
* Merge Assemblies of macaque and ebov Virus to generate one merged assembly.
*/
process merge_assemblies {

    storeDir "${params.output_dir_preliminary}/reference_assembly"

    input:
    file rheMac from rhesus_genome_channel
    file ebov from ebov_genome_channel

    output:
    set file("${params.prefix}.fa"), file("${params.prefix}.fa.fai") into (merged_assembly,fasta_reference_channel,  reference_assembly_channel,  merged_assembly_for_dictionary )

    script:
    """
    cat ${rheMac} > ${params.prefix}.fa
    sed 's/KU182905.1/EBOV_Kikwit/g' ${ebov} >> ${params.prefix}.fa
    samtools faidx ${params.prefix}.fa
    """

}

/*
* Concatenate ensembl release annotation with EBOV gene annotation
*/
process merge_annotations{
  storeDir "${params.output_dir_preliminary}/gene_annotations"

  input:
  file ebov_annotation from ebov_annotation_channel
  file(rheMac_annotation) from rheMac_annotation_channel

  output:
  file("${params.prefix}.gtf") into (merged_annotation_channel,  merged_annotation_ch, merged_annotation_ch_2  )

  script:
  """
  cat  ${rheMac_annotation} ${ebov_annotation} > ${params.prefix}.gtf
  """

}

//need to add the other part
process convert_gtf_to_bed12{
  storeDir "${params.output_dir_preliminary}/gene_annotations"

  input:
  file merged_annotation from merged_annotation_ch
  file gtfToGenePred_script from gtfToGenePred_script_ch
  file genePredToBed_script from genePredToBed_script_ch

  output:
  set file("${params.prefix}.bed"), file("${params.prefix}.bed12") into bed_channel

  script:
  """
  ./${gtfToGenePred_script} ${params.prefix}.gtf ${params.prefix}.bed
  ./${genePredToBed_script} ${params.prefix}.bed ${params.prefix}.bed12
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

  storeDir "${params.output_dir_preliminary}/indexes/hisat2"
  cpus 2

  input:
  set file(assembly), file(fai) from merged_assembly

  output:
  file "${params.prefix}.*" into hisat2_indexes

  script:
  """
  hisat2-build-s ${assembly} ${params.prefix} -p ${task.cpus}
  """
}


/*
* Create dictionary for merged assembly with picard tools
*/
process create_dictionary{

  storeDir "${params.output_dir_preliminary}/reference_assembly/"

  input:
  set file(assembly), file(fai) from merged_assembly_for_dictionary

  output:
  file "${params.prefix}.dict" into dictionary_channel

  script:
  """
  picard-tools CreateSequenceDictionary R=${assembly} O=${params.prefix}.dict
  """
}


/*
* Extract exons and splice sites for annotated genes
*/
process extract_exons_ss{
 storeDir "${params.output_dir_preliminary}/gene_annotations"

 input:
 file merged_annotation from merged_annotation_channel
 file extract_ss from extract_ss_script
 file extract_exon from extract_exon_script
 output:
 set file("${params.prefix}.exon.txt"), file("${params.prefix}.ss.txt") into extracted_exons_ss_channel

 shell:
 '''
 ./!{extract_ss} !{merged_annotation} > !{params.prefix}.ss.txt
 ./!{extract_exon} !{merged_annotation} > !{params.prefix}.exon.txt
 '''
}


log.info "=============================================="
log.info " RNA-Seq Analysis  "
log.info "=============================================="



// Channel duplications
merged_annotation_ch_2.into{ gtfChannel1; gtfChannel2; gtfChannel3; gtfChannel4; gtfChannel5; gtfChannel6; ref_gtf_channel; ref_gtf_channel_2}



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
      file(fastq_pair) from fastq_files_for_qc

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
* The indexes were computed WITHOUT the help of known splice sites
* the known splice sites are now used as input for the mapper.
*/
process mapping_hisat{

  label 'big_mem'
  cpus 1
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

  input:
  file indexes from hisat2_indexes.collect()
  set file(exon), file(ss) from extracted_exons_ss_channel.collect()
  set dataset_name, tissue, dayPostInfection, sample,lane, complete_id,
      file(fastq_pair)  from fastq_files_for_mapping

  output:
  set lane, dataset_name, dayPostInfection, tissue, sample, complete_id,
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
  hisat2 -p ${task.cpus} -x ${params.prefix} -1   ${fastq_pair[0]}  \
                      -2  ${fastq_pair[1]}  --known-splicesite-infile ${ss} \
                      --novel-splicesite-outfile ${complete_id}.novel_ss.txt \
                      --downstream-transcriptome-assembly \
                      --time --summary-file ${complete_id}.hisat2_summary.txt \
                      --rna-strandness ${params.strandness} > ${complete_id}.sam
  """

}



/*
* STEP 3.2: Sort and index the mapped bam files
*
* The sorting is done with picard as we encountered some problems when running
* samtools and afterwards other picardTools on top of the sorted bam.
*/
process sort_bam{

  cpus 12
  tag "${complete_id}"
  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"


  input:
  set lane, dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(novel_ss),file(summary),
      file(sam) from mapped_sam

  output:
  set complete_id,lane, dataset_name, dayPostInfection, tissue, sample,
      file(summary),
      file(sam), file("${complete_id}.bam")   into (sorted_and_index_bam, sorted_indexed_bams_for_stats, hisat2_bams)
  file("${complete_id}.bam") into sorted_and_index_bam_2

  script:
  """
  samtools sort ${sam} -T $TMPDIR/${complete_id}.tmp -@ ${task.cpus} -O bam -o ${complete_id}.bam
  """
}

/*
* STEP 4: Merge bam alignments -  Add UMIs
*
* After having converted the unmapped bam files into fastq we have lost
* relevant information like UMIs info.
* We now merge the unmapped and the mapped bam files to retrieve information
* about the UMIs.
*
*  A command-line tool for merging BAM/SAM alignment info from a third-party
*  aligner with the data in an unmapped BAM file, producing a third BAM file that
*  has alignment data (from the aligner) and all the remaining data from the
*  unmapped BAM. Quick note: this is not a tool for taking multiple sam files
* and creating a bigger file by merging them.
*
*/

// We join the mapped and unmapped channel by complete_id.
// See how channel were built at the beginning - complete id is always in the -2 position

if( "${params.umi}" == "true" ){
  unmapped_and_mapped_bams = sorted_and_index_bam.combine(unmapped_bams, by:0)

  process add_unmapped_bam{

    tag "${complete_id}"
    storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"
    label 'big_mem'

    input:
    set complete_id,lane, dataset_name, dayPostInfection, tissue, sample,
        file(summary),
        file(sam), file(bam),
        file(unmapped_bam) from unmapped_and_mapped_bams
    set file(assembly), file(fai) from reference_assembly_channel.collect()
    file dict from dictionary_channel.collect()

    output:
    file ("${complete_id}.UMI.bam") into merged_bam_alignments_channel

    script:
    """
    picard-tools  MergeBamAlignment UNMAPPED=${unmapped_bam} ALIGNED=${bam} \
                  O=${complete_id}.UMI.bam R=${assembly} \
                  SO=coordinate ALIGNER_PROPER_PAIR_FLAGS=true MAX_GAPS=-1 \
                  ORIENTATIONS=${params.strandness} CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT \
                  PAIRED_RUN=true
    """
  }
}
else{
  sorted_and_index_bam_2.set{merged_bam_alignments_channel}
}


// // In this channel we group the files by the complete_id after removing the lane.
// // We have now a set of bams per sample from all the lanes.
merged_bam_alignments_channel.map { file ->
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


/*
* Step 5 Merge lanes
* Each sample has 4 lanes.
* We need to merge all of them into one.
*
* Aaron's mail
* We prepared both the old and new Zyagen samples at the same time with the s
* same library prep. We introduce UMIs at the adaptor ligation step, prior to any PCR.
* The bam files in the recent run folder represent a single library prep per
* sample, sequenced over four lanes, for both the old and new Zyagen samples.
* So we should merge bams for all four lanes, and then
*/
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
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(merged_bam), file(bai_mereged_bam) from merged_bylanes

  output:
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file("${complete_id}.UMI.f3.q60.bam"), file("${complete_id}.UMI.f3.q60.bam.bai") into filtered_bams
  script:
  """
  samtools view -b -f3 -q 60  ${merged_bam} > ${complete_id}.UMI.f3.q60.bam
  samtools index ${complete_id}.UMI.f3.q60.bam
  """
}

filtered_bams_ch = Channel.create()
/*
* STEP 7B : Dedup duplicate using UMIs
*  --paired BAM is paired end - output both read pairs. This will also force the use of
*  the template length to determine reads with the same mapping coordinates.
* -- directional (default)
*  Identify clusters of connected UMIs (based on hamming distance threshold)
*  and umi A counts >= (2* umi B counts) - 1. Each network is a read group.
*/
if( "${params.umi}" == "true" ){
  process dedupUmi{

    tag "${complete_id}"
    label 'big_mem'
    cpus 48
    storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample"

    input:
    set dataset_name, dayPostInfection, tissue, sample, complete_id, file(filtered_bam), file(filtered_bai) from filtered_bams

    output:
    set dataset_name, dayPostInfection, tissue, sample, complete_id,
        file("${complete_id}.UMI.f3.q60.umi_dedup.bam"),
        file("${complete_id}.UMI.f3.q60.umi_dedup.bam.bai") into filtered_bams_ch

    script:
    """
    mkdir umi_logs
    umi_tools dedup --stdin=${filtered_bam} --log=${complete_id}.umi_tools.dedup.log --output-stats=${complete_id} --extract-umi-method=tag --umi-tag=RX --paired --method directional  > ${complete_id}.UMI.f3.q60.umi_dedup.bam
    mv ${complete_id}.umi_tools.dedup.log umi_logs/
    mv ${complete_id}_per* umi_logs/
    mv ${complete_id}_edit* umi_logs/
    samtools index ${complete_id}.UMI.f3.q60.umi_dedup.bam
    """

  }
}else{

  filtered_bams.into(filtered_bams_ch)
}

filtered_bams_ch.into{ filtered_merged_bams_1;  filtered_merged_bams_2;  filtered_merged_bams_3;   filtered_merged_bams_4;   filtered_merged_bams_5}


process runRSeQC{

  storeDir "${params.output_dir}/03_hisat/$dataset_name/$tissue/$dayPostInfection/$sample/rseqc"

  input:
  set  dataset_name, dayPostInfection, tissue, sample, complete_id,
      file(bam),
      file(bai) from filtered_merged_bams_2
  set file(bed), file(bed12) from  bed_channel.collect()

  output:
  file "*" into read_distribution_channel

  script:
  """
  read_distribution.py -i ${bam} -r ${bed12} > ${complete_id}.UMI.f3.q60.read_distribution.txt
  infer_experiment.py -i ${bam} -r ${bed12} > ${complete_id}.infer_experiment.txt
  junction_annotation.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12}
  bam_stat.py -i ${bam} > ${complete_id}.bam_stat.txt
  junction_saturation.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12} > ${complete_id}.junction_annotation_log.txt
  inner_distance.py -i ${bam} -o ${complete_id}.rseqc -r ${bed12}
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
  set dataset_name, dayPostInfection, tissue, sample, complete_id, file(bam), file(bai) from filtered_merged_bams_1
  file gtf from gtfChannel4.collect()

  output:
  val dataset_name into dataset_name_ch
  file "*.stringtie.gtf" into stringTie_channel

  script:
  file_prefix = get_file_name_no_extension(bam.name)
  """
  stringtie ${bam} -G ${gtf} -o ${file_prefix}.stringtie.gtf --fr -p ${task.cpus}
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
  file(stringtie_gtfs) from stringTie_channel.collect()
  file reference_gtf from gtfChannel5

  output:
  file "stringtie_merged_reference_guided.gtf" into (merged_denovo_assmebly, merged_de_novo_assembly_2)

  script:
  """
  stringtie --merge -p ${task.cpus} -F 1.0 -o stringtie_merged_reference_guided.gtf -G ${reference_gtf} ${stringtie_gtfs}
  """
}


/*
* Obtain stats over StringTie output
*/
process gffCompare2{

  storeDir "${params.output_dir}/04b_gffCompare"

  input:
  file merged_gtf from merged_denovo_assmebly
  file reference_gtf from gtfChannel6
// val dataset_name from dataset_name_ch2

  output:
  file("merged*") into gff_compare_output_channel2

  script:
  """
  gffcompare -R -r ${reference_gtf} -o merged ${merged_gtf}
  """
}


merged_de_novo_assembly_2.into{ merged_de_novo_assembly_3; merged_de_novo_assembly_4; merged_de_novo_assembly_5}



log.info "=============================================="
log.info "            lncRNA annotation"
log.info "=============================================="

process feelnc_filter{

  storeDir "${params.output_dir}/05_lncrnaAnnotation_monoexonic_included/feelnc_gencode_linc"

  input:
  file merged_gtf from merged_de_novo_assembly_3
  file reference_gtf from ref_gtf_channel

  output:
  file("candidate_lncRNA.gtf") into candidates

  script:
  """
  FEELnc_filter.pl -i ${merged_gtf} \
                   -a ${reference_gtf} \
                   --monoex=-1 \
                   -b transcript_biotype=protein_coding > candidate_lncRNA.gtf
  """
}

process feelnc_codpot{
  cpus 48
  storeDir "${params.output_dir}/05_lncrnaAnnotation_monoexonic_included/feelnc_gencode_linc"

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

  storeDir "${params.output_dir}/05_lncrnaAnnotation_monoexonic_included/feelnc_gencode_linc/01_gffcompare"

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

// ---------------------------------------
// Explorative - not real part of the pipeline
// ----------------------------------------


// process getHTseqCountMerged{
//
//   cpus 8
//   storeDir "${params.output_dir}/06_counts/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"
//
//   input:
//   file gtf from merged_de_novo_assembly_4.collect()
//   set dataset_name, dayPostInfection, tissue, sample, complete_id,
//       file(bam),
//       file(bai) from filtered_merged_bams_3
//
//   output:
//   file "${bam_prefix}.HTseq.gene_counts.tab" into htseqCountsMerged_channel
//
//   script:
//   bam_prefix = get_file_name_no_extension(bam.name)
//   """
//   # Calc the counts for the umi_dedup
//   htseq-count -f bam -r name -s ${params.htseqsense} -t exon -i gene_id ${bam} ${gtf} > ${bam_prefix}.HTseq.gene_counts.tab
//   """
// }
//
//
// process getHTseqCountRef{
//
//   cpus 8
//   storeDir "${params.output_dir}/06_counts_ref/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"
//
//   input:
//   file gtf from gtfChannel1.collect()
//   set dataset_name, dayPostInfection, tissue, sample, complete_id,
//       file(bam),
//       file(bai) from filtered_merged_bams_4
//
//   output:
//   file "${bam_prefix}.HTseq.gene_counts.tab" into htseqCountsRef_channel
//
//   script:
//   bam_prefix = get_file_name_no_extension(bam.name)
//   """
//   # Calc the counts for the umi_dedup
//   htseq-count -f bam -r name -s ${params.htseqsense} -t exon -i gene_id ${bam} ${gtf} > ${bam_prefix}.HTseq.gene_counts.tab
//   """
// }
//
// // I merge the novel with the reference
// params.gtf_ref_merged =  "/gpfs/projects/bsc83/Data/Ebola/01_Ebola-RNASeq_all/03_novel_lncrnas/00_gtf_external/rheMac10_EBOV-Kikwit_and_novel.gtf"
// Channel.fromPath("${params.gtf_ref_merged}").set{ gtfANDnovel }
//
// process getHTseqCountALL{
//
//   cpus 24
//   storeDir "${params.output_dir}/06_counts_all/$dataset_name/$tissue/$dayPostInfection/$sample/htseq_counts"
//
//   input:
//   file gtf from gtfANDnovel.collect()
//   set dataset_name, dayPostInfection, tissue, sample, complete_id,
//       file(bam),
//       file(bai) from filtered_merged_bams_5
//
//   output:
//   file "${bam_prefix}.HTseq.gene_counts.tab" into htseqCountsALL_channel
//
//   script:
//   bam_prefix = get_file_name_no_extension(bam.name)
//   """
//   # Calc the counts for the umi_dedup
//   htseq-count -f bam -r name -s ${params.htseqsense} -t exon -i gene_id ${bam} ${gtf} > ${bam_prefix}.HTseq.gene_counts.tab
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
