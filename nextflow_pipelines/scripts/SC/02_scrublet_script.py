import scrublet as scr
import scipy.io
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('input')
parser.add_argument('name_mask')
parser.add_argument('name_score')


args = parser.parse_args()
#counts_matrix = scipy.io.mmread("/home/luisas/Desktop/cluster/data/RObjects/sparse_invivo.mtx").T.tocsc()
counts_matrix = scipy.io.mmread(args.input).T.tocsc()
#print('Counts matrix shape: {} rows, {} columns'.format(counts_matrix.shape[0], counts_matrix.shape[1]))
scrub = scr.Scrublet(counts_matrix, expected_doublet_rate=0.06)
doublet_scores, predicted_doublets = scrub.scrub_doublets(min_counts=2, 
                                                          min_cells=3, 
                                                          min_gene_variability_pctl=85, 
                                                          n_prin_comps=30)
scrub.call_doublets(threshold=0.25)
np.savetxt(args.name_mask, scrub.predicted_doublets_, fmt='%s')
np.savetxt(args.name_score, scrub.predicted_doublets_, fmt='%s')
