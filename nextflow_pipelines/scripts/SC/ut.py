import os
import pandas as pd
import numpy as np
import scanpy as sc


###########################################################################################
### Functions for loading files from google cloud storage, and H5AD files in particular ###

from tempfile import NamedTemporaryFile
import subprocess

def save_adata(adata, filepath, ext='.h5ad', gcs=False):
    if gcs:
        temp = NamedTemporaryFile(suffix=ext, delete=False)
        temp.close()
        sc.write(temp.name, adata)
        subprocess.call('gsutil -m cp %s %s' % (temp.name, filepath), shell=True)
        subprocess.call('rm %s' % temp.name, shell=True)
    else:
        sc.write(filepath, adata)
    
    
    
def read_adata(filepath, ext='.h5ad', gcs=False):
    if gcs:
        temp = NamedTemporaryFile(suffix=ext, delete=False)
        temp.close()
        subprocess.call('gsutil -m cp %s %s' % (filepath, temp.name), shell=True)
        adata = sc.read(temp.name)
        subprocess.call('rm %s' % temp.name, shell=True)
    else:
        adata = sc.read(filepath)
    return(adata)

#
