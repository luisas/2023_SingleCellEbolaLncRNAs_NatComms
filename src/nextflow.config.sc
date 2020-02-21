process{

  container = 'file:///gpfs/projects/bsc83/containers/singularity/dropseq.simg'

  withName:FastqToBam{
    container = 'file:///gpfs/projects/bsc83/containers/singularity/rnaseqnew.simg'
  }
  withName:create_star_indexes{
    container = 'file:///gpfs/projects/bsc83/containers/singularity/rnaseq.simg'
  }
  withName:STAR{
    container = 'file:///gpfs/projects/bsc83/containers/singularity/rnaseq.simg'
  }
  withLabel:rnaseq{
    container = 'file:///gpfs/projects/bsc83/containers/singularity/rnaseqnew.simg'
  }
  withName:CellSelection{
    container = 'file:///gpfs/projects/bsc83/containers/singularity/rdropbead.simg'
  }
  withName:runRSeQC{
   container = 'file:///gpfs/projects/bsc83/containers/singularity/rseqc.simg'
  }


}

singularity {
    enabled = true
}
