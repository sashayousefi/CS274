import sys
import numpy as np
import pandas as pd
import chemoUtils
import itertools
import pvalue
import networkx as nx
from collections import deque
import matplotlib.pyplot as plt

class Network(object):
    """
    Network class to attain a graphical representation of our protein nodes. A network where 
    the nodes are proteins and edges are drawn between pairs for which the bootstrap 
    p-value is less than or equal to 0.05
    """
    def __init__(self, drug_file, target_file, protein_nodes):
        """
            Function to initialize objects holding the necessary drug/protein info for our 
            graphical representation
        """
        self.drug_file = drug_file
        self.target_file = target_file
        #read in protein_nodes into a dataframe with protein name/id/indication
        self.protein_nodes = pd.read_csv(protein_nodes, sep = ',')
        #list of proteins in protein nodes
        self.protein_lst = self.protein_nodes.loc[:, 'uniprot_accession'].values
        #initialize a p_value dictionary with lists of proteins in significant pairs for our network_edgelist
        self.pval_df = {'prot1':[], 'prot2':[]}

    def generate_network(self):
        """
        Generate a graphical representation of our ligands 
        --------
        Params:
            - None: use init defined protein/drug values
        Updates:
            - self.pval_df: dataframe with 
        """
        drug_pairs = itertools.combinations(self.protein_lst, 2)
        for pairs in drug_pairs:
            drug_i = pairs[0]
            drug_j = pairs[1]
            pval_obj = pvalue.pvalue(n= 500, r = 214, drug_file = self.drug_file, target_file = self.target_file, proteinA = drug_i, proteinB = drug_j)
            pval = pval_obj.bootstrap(print_pval = False)
            if pval <= 0.05:
                self.pval_df['prot1'].append(drug_i)
                self.pval_df['prot2'].append(drug_j)
        self.pval_df = pd.DataFrame.from_dict(self.pval_df)

    def generate_edgelist(self):
        self.generate_network()
        self.pval_df.to_csv('network_edgelist.txt', header=False, index = False, sep = ' ')


def main():
    # input variables
    drug_file = sys.argv[1]
    target_file = sys.argv[2]
    protein_nodes = sys.argv[3]

    network_ex = Network(drug_file, target_file, protein_nodes)
    network_ex.generate_edgelist()


if __name__== "__main__":
    main()
