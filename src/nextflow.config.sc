process{

  container = 'file:///gpfs/projects/bsc83/containers/singularity/dropseq.simg'

  withName:CellSelection{
    container = 'file:///gpfs/projects/bsc83/containers/singularity/rdropbead.simg'
  }


}

singularity {
    enabled = true
}
