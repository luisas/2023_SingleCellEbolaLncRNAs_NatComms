process{

  container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/dropseq.simg'

  withName:Sam2Fastq{
    container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/rnaseqnew.simg'
  }
  withName:create_star_indexes{
    container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/rnaseq.simg'
  }
  withName:STAR{
    container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/rnaseq.simg'
  }
  withName:STAR_unpaired{
    container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/rnaseq.simg'
  }
  withLabel:rnaseq{
    container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/rnaseqnew.simg'
  }
  withName:CellSelection{
    container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/rdropbead.simg'
  }
  withName:runRSeQC{
   container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/rseqc.simg'
  }
  withName:FPKM{
   container = 'file:///gpfs/projects/bsc83/utils/containers/singularity/rseqc.simg'
  }


}

singularity {
    enabled = true
}
