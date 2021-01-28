
log.info "=============================================="
log.info "            Benchmark: LncRNA annotation      "
log.info "=============================================="

/*
 * This pipeline is built to benchmark the performance of annotation tools
 * such as CPC2, CPAT and CNIT, on the reference annotation of Macaque to
 * assess their performance.
 */

// CONFIG FILE: ../configs/nextflow.annotation.config

// Input directory
params.input_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/"

// Output directory
output_dir_sub="${params.input_dir}/03_novel_lncRNAs_list/99_benchmark_annotation/"

// Reference genome ( Macaca Mulatta rheMac10 )
reference_assembly = Channel.fromPath("/gpfs/projects/bsc83/Data/assemblies/ensembl/release-100/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa").collect()

// Reference annotation ( Macaca Mulatta rheMac10)
gtf = Channel.fromPath("/gpfs/projects/bsc83/Data/gene_annotation/ensembl_release100/rheMac10/Macaca_mulatta.Mmul_10.100.gtf")

// Script for conversion of gtf to bed format
gtf2bed = Channel.fromPath("${baseDir}/scripts/gtf2bed0.R").collect()

// Extra (default) files needed by CPAT for obtaining the prediction
cpatmodel= Channel.fromPath("${params.input_dir}/01_PreliminaryFiles_rheMac10/Human_logitModel.RData").collect()
cpathexamer= Channel.fromPath("${params.input_dir}/01_PreliminaryFiles_rheMac10/Human_Hexamer.tsv").collect()
cpat_files = Channel.fromPath("${params.input_dir}/01_PreliminaryFiles_rheMac10/cpat*").collect()



/*
 * 1) Convert GTF to BED
 */

process gtf2bed{
   storeDir "${output_dir_sub}/00_preliminaryfiles/"

   input:
   file gtf2bed
   file candidates from gtf

   output:
   file("${candidates.baseName}.bed12") into candidatesbed

   script:
   """
   Rscript ${gtf2bed} ${candidates} ${candidates.baseName}.bed12
   """
}

/*
 * 2) Obtain FASTA file
 */
process getFasta{
   storeDir "${output_dir_sub}/00_preliminaryfiles/"

   input:
   file reference_assembly
   file candidates from candidatesbed

   output:
   file("${candidates.baseName}.fa") into candidates

   script:
   """
   bedtools getfasta -fi ${reference_assembly} -bed ${candidates} -name -fo ${candidates.baseName}.fa -s
   """
}

candidates.into{candidatesfa1; candidatesfa2; candidatesfa3; }


/*
 * 3) Run prediction with CPC2
 */
process CPC2{
  storeDir "${output_dir_sub}/01_predictions/CPC2"
   input:
  file candidates from candidatesfa1
   output:
   file("*") into cpc2_all
   file("cpc2_pred.txt") into cpc2

   script:
   """
   python2 \${CPC_HOME}/CPC2.py -i ${candidates} -o cpc2_pred
   """
}

/*
 * 4) Run prediction with CPAT
 */
process CPAT{
   storeDir "${output_dir_sub}/01_predictions/CPAT"

   input:
   file candidates from candidatesfa2
   file cpatmodel
   file cpathexamer
   file cpat_files

   output:
   file("*") into cpat_all
   file("cpat_pred.ORF_prob.tsv") into cpat

   script:
   """
   cpat.py -g ${candidates} -o cpat_pred  -d ${cpatmodel} -x ${cpathexamer}
   """
}


/*
 * 5) Run prediction with CNIT
 */
process CNIT{
    storeDir "${output_dir_sub}/01_predictions/CNIT"

    input:
    file candidates from candidatesfa3

    output:
    file("*") into cnit_all
    file("cnit_pred/CNCI2.index") into cnit
    script:
    """
    python \${CNIT_HOME}/CNCI2.py -f ${candidates} -o cnit_pred -m "ve"
    """
}
