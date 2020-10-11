
// Filter is performed on the correct RheMac gtf

params.output_dir = "/gpfs/projects/bsc83/Data/Ebola/01_bulk_RNA-Seq_lncRNAs_annotation/"
params.dataset= "02_RNA-Seq_ribodepl"
params.prefix_label = "ribodepleted"
output_dir_sub="03_novel_lncRNAs_list/99_benchmark_annotation/"

candidates = Channel.fromPath("/gpfs/projects/bsc83/Data/assemblies/ensembl/release-100/rheMac10/Macaca_mulatta.Mmul_10.dna.toplevel.fa").collect()

candidates.into{candidatesfa1; candidatesfa2; candidatesfa3; }

cpatmodel= Channel.fromPath("${params.output_dir}/01_PreliminaryFiles_rheMac10/Human_logitModel.RData").collect()
cpathexamer= Channel.fromPath("${params.output_dir}/01_PreliminaryFiles_rheMac10/Human_Hexamer.tsv").collect()
cpat_files = Channel.fromPath("${params.output_dir}/01_PreliminaryFiles_rheMac10/cpat*").collect()


log.info "=============================================="
log.info "            LncRNA annotation"
log.info "=============================================="

process CPC2{
  storeDir "${params.output_dir}/${output_dir_sub}/01_predictions/CPC2"
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


process CPAT{
   storeDir "${params.output_dir}/${output_dir_sub}/01_predictions/CPAT"

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



process CNIT{
    storeDir "${params.output_dir}/${output_dir_sub}/01_predictions/CNIT"

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
