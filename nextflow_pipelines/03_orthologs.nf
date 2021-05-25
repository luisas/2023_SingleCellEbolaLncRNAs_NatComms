

//  To run locally

//  Also Download Slnky
// TODO:
// wget https://github.com/slncky/slncky/archive/refs/heads/master.zip
// unzip master.zip
// rm master.zip

// Download dependencies
// From  https://www.dropbox.com/s/jifnrr3kjq95w64/annotations.tar.gz?dl=0 into slinky master dir
// tar -xf annotations.tar.gz
// rm annotations.tar.gz

// Also TODO
// cp /home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/00_preliminaryfiles/rheMac10.fa ./slncky-master/annotations/
// cp /home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/00_preliminaryfiles/rheMac10.fa.fai ./slncky-master/annotations/
// cp /home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/00_preliminaryfiles/rheMac10ToHg38.over.chain.gz ./slncky-master/annotations/
// cp /home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/00_preliminaryfiles/annotations.config ./slncky-master/

params.dirData = "/home/luisas/Desktop/cluster/data/"
params.slnkydir = "/home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/slncky-master/"
Channel.fromPath("${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/rheMac10_EBOV_and_novel_genenames.gtf")
       .into{ gtfChannel; GtfChannel2; GtfChannel3;   }

prep_subfiles = Channel.fromPath("${baseDir}/scripts/createSubfiles.R").collect()
convert =  Channel.fromPath("${baseDir}/scripts/convert_annotation.R").collect()
gtfToGenePred = Channel.fromPath("${baseDir}/scripts/gtfToGenePred").collect()
genePredToBed = Channel.fromPath("${baseDir}/scripts/genePredToBed").collect()
config = Channel.fromPath("${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/00_preliminaryfiles/annotations.config").collect()
dirslnky = Channel.fromPath("${params.slnkydir}")
sourcefile = Channel.fromPath("${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/00_preliminaryfiles/ensemblSource.txt").collect()
chrom = Channel.fromPath("${params.slnkydir}/annotations/chromosomes_in_fasta.tsv").collect()



params.prefix = "rheMac10"

params.out_preliminary = "${params.dirData}/01_bulk_RNA-Seq_lncRNAs_annotation/03_novel_lncRNAs_list/"
// 1 . Create other subfiles (Separate gtfs for separate genes biotypes )
// gtf <- import(args[1])
// sourceTable <- read.table(args[2])
// outdir <- args[3]
// prefix <- args[4]
// ref_novel_full<- import(args[5])
//Rscript /home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/slncky-master/annotations/createSubfiles.R

process covertAnnotationUCSC{

  storeDir "${params.out_preliminary}"

  input:
  file convert
  file gtf from GtfChannel2

  output:
  file("rheMac10_EBOV_and_novel_genenames_UCSC.gtf") into (refFullChannel,refFullChannel2)

  script:
  """
  Rscript ${convert} ${gtf} "UCSC" "${params.out_preliminary}/rheMac10_EBOV_and_novel_genenames_UCSC.gtf"
  """

}


process Prepare_Subfiles{
  publishDir "${params.slnkydir}/annotations", mode: 'copy', overwrite: false
  input:
  file prep_subfiles
  file gtf from gtfChannel
  file reffull from refFullChannel
  file gtfToGenePred
  file genePredToBed
  file sourcefile
  file chrom


  output:
  file "*" into prep_subfiles_channel
  val "done" into prepdone
  file "rheMac10_EBOV_and_novel_genenames_UCSC.contigsfiltered.gtf" into contigFull

  script:
  """
  Rscript ${prep_subfiles} ${gtf} ${sourcefile} "." ${params.prefix} ${reffull} ${chrom}
  ./gtfToGenePred ${params.prefix}.ensGene.lnc.gtf ${params.prefix}.ensGene.lnc.genepred
  ./gtfToGenePred ${params.prefix}.ensGene.lnc.gtf ${params.prefix}.ensGene.lnc.genepred
  ./genePredToBed ${params.prefix}.ensGene.lnc.genepred ${params.prefix}.ensGene.lnc.bed

  ./gtfToGenePred ${params.prefix}.ensGene.mirna.gtf ${params.prefix}.ensGene.mirna.genepred
  ./genePredToBed ${params.prefix}.ensGene.mirna.genepred ${params.prefix}.ensGene.mirna.bed

  ./gtfToGenePred ${params.prefix}.ensGene.PC.gtf ${params.prefix}.ensGene.PC.genepred
  ./genePredToBed ${params.prefix}.ensGene.PC.genepred ${params.prefix}.ensGene.PC.bed

  ./gtfToGenePred ${params.prefix}.ensGene.Pseudogenes.gtf ${params.prefix}.ensGene.Pseudogenes.genepred
  ./genePredToBed ${params.prefix}.ensGene.Pseudogenes.genepred ${params.prefix}.ensGene.Pseudogenes.bed

  ./gtfToGenePred ${params.prefix}.ensGene.PC.gtf ${params.prefix}.ensGene.PC.genepred
  ./genePredToBed ${params.prefix}.ensGene.PC.genepred ${params.prefix}.ensGene.PC.bed

  ./gtfToGenePred ${params.prefix}.ensGene.snorna.gtf ${params.prefix}.ensGene.snorna.genepred
  ./genePredToBed ${params.prefix}.ensGene.snorna.genepred ${params.prefix}.ensGene.snorna.bed
  """
}

process convert{

  storeDir "${params.out_preliminary}"

  input:
  file gtfToGenePred
  file genePredToBed
  file gtf from GtfChannel3
  file reffull from contigFull


  output:
  file "*" into conversions
  file "${reffull.baseName}.bed" into ref_contigsfiltered_bed

  script:
  """
  ./gtfToGenePred ${gtf} ${gtf.baseName}.genepred
  ./genePredToBed ${gtf.baseName}.genepred ${gtf.baseName}.bed

  ./gtfToGenePred ${reffull} ${reffull.baseName}.genepred
  ./genePredToBed ${reffull.baseName}.genepred ${reffull.baseName}.bed
  """

}


// cut -f1 rheMac10.fa.fai > chromosomes_in_fasta.tsv
//./slncky.v1.0 -2 /home/luisas/Desktop/cluster/data/01_bulk_RNA-Seq_lncRNAs_annotation/05_orthologs/slncky-master/rheMac10_EBOV_and_novel_genenames_UCSC.contigsfiltered.bed mmul10 out




// For the moment leave out of pipeline and run manually - on a hurry :S
// process slnky{
//
//   publishDir "${params.slnkydir}", mode: 'copy', overwrite: true
//
//   input:
//   file dirslnky
//   file config
//   file bed from ref_contigsfiltered_bed
//   val done from prepdone
//   file annot
//   file chain
//
//   output:
//   file "*" into slnkyoutput
//
//   script:
//   """
//   cd ${dirslnky}
//   ./slncky.v1.0 -c ${config} -2  ${bed} mmul10 out
//   """
// }
