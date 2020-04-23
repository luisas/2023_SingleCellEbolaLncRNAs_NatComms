// BaseFolders
params.dirData = "/gpfs/projects/bsc83/Data/Ebola"
params.dirProj = "/gpfs/projects/bsc83/Projects/Ebola"


//params.dirData = "/Users/luisasantus/Desktop/cluster/data"
//params.dirProj = "/Users/luisasantus/Desktop/cluster/proj"
params.output_dir_preliminary = "${params.dirData}/01_Ebola-RNASeq_all/01_PreliminaryFiles_rheMac10/"

params.output_dir = "${params.dirData}/01_Ebola-RNASeq_all/03_scRNA-Seq_complete"



Channel.fromPath("${params.output_dir}/03_star/*/*/*/*/*/*.sam")
              .ifEmpty('bam files directory is empty')
              .map{tuple(it.toString().split('/')[8,9,10,11,12,13].join("/"),
              it.toString().split('/')[-1].split("\\.").init().join("."),
                         it)}
              .into{ bams_1; bams_2; bams_3;  }

params.scripts="${baseDir}/scripts/"
gtfToGenePred_script_ch = Channel
                      .fromPath("${params.scripts}/gtfToGenePred")
genePredToBed_script_ch = Channel
                      .fromPath("${params.scripts}/genePredToBed")

Channel.fromPath("${params.dirData}/01_Ebola-RNASeq_all/03_novel_lncrnas/02_final_catalogue/rheMac10_EBOV-Kikivit_and_novelcatalogue_with_names.gtf")
       .into{ GtfChannel;GtfChannel2; GtfChannel3;GtfChannel4; GtfChannel5;   }
//Channel.fromPath("${params.dirData}/01_Ebola-RNASeq_all/03_novel_lncrnas/02_final_catalogue/03_liftover/rheMac8_EBOV-Kikivit_and_novelcatalogue_with_names_ensembl.gtf")
//       .into{ GtfChannel;GtfChannel2; GtfChannel3;GtfChannel4; GtfChannel5;   }

Channel.fromPath("${params.output_dir_preliminary}/gene_annotations/rheMac10_EBOV-Kikwit.gtf")
        .into{ GtfNotNovel;   }

Channel.fromPath("${params.output_dir_preliminary}/reference_assembly/rheMac10_EBOV-Kikwit.fa")
       .into{ ReferenceChannel;ReferenceChannel2;   }

Channel.fromPath("${params.dirData}/01_Ebola-RNASeq_all/01_PreliminaryFiles_rheMac10/reference_assembly/rheMac10_EBOV-Kikwit.dict")
       .into{ DictChannel; DictChannel2;DictChannel3; DictChannel4;  }

params.strandness = "FR"



cell_selection_script = Channel.fromPath("${baseDir}/scripts/cellselection.R").collect()
// We need to estimate how manty cells we want to extract

bams_2.subscribe{ println "$it" }

/*
process compress{

  cpus 12
  storeDir "${params.output_dir}/$path"

  input:
  set path, filename, file(sam) from bams_1

  output:
  set path, filename, file("${filename}.bam") into bamba

  script:
  """
  samtools view -Sb  $sam  >  ${filename}.bam
  """
}
*/
