import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2
from distributed import LocalCluster, Client
import sys


if __name__ == '__main__':

    # ex_matrix is a DataFrame with gene names as column names
    #ex_matrix = pd.read_csv("/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/matrix_final_myeloids_filtered.csv", sep=',')
    ex_matrix = pd.read_csv(sys.argv[1], sep=',')
    print("1) Gene expression matrix has been succesfully read")

    network = grnboost2(expression_data=ex_matrix)

    print("2) Network has been succesfully created")

    network.to_csv(sys.argv[2], sep='\t', index=False, header=False)
    #network.to_csv('/gpfs/projects/bsc83/Data/Ebola/02_scRNA-Seq_PBMCs/01_scRNA-Seq_inVivo_rhemac10/05_RObjects/09_graph/output_myeloids_filtered.tsv', sep='\t', index=False, header=False)
    print("3) DONE!")
