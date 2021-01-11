import pandas as pd
from arboreto.utils import load_tf_names
from arboreto.algo import grnboost2



if __name__ == '__main__':
    # ex_matrix is a DataFrame with gene names as column names
    ex_matrix = pd.read_csv("/gpfs/projects/bsc83/Data/Ebola/matrix.csv", sep=',')


    network = grnboost2(expression_data=ex_matrix)

    network.to_csv('output.tsv', sep='\t', index=False, header=False)

