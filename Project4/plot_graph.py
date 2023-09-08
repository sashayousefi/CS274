import sys
import numpy as np
import pandas as pd
import chemoUtils
import itertools
import pvalue
import networkx as nx
from collections import deque
import matplotlib.pyplot as plt

class NetworkGraph(object):
    """
    NetworkGraph class to draw a graphical representation of our protein nodes.
    """
    def __init__(self, network_edgelist, protein_nodes, output_file):
        """
            Function to initialize objects holding the necessary drug/protein info for our 
            graphical representation
        """
        self.network_edgelist = network_edgelist
        self.protein_nodes = pd.read_csv(protein_nodes, sep = ',')
        self.output_file = output_file
        #mapping of indication to color
        self.color_mapping = {"bp":"red", "bp;cholesterol":"green", "bp;cholesterol;diabetes":"blue", "bp;diabetes":"purple"}
        #mapping of uniprot_accession to uniprot_id
        self.names_mapping = {}
        for i in np.arange(len(self.protein_nodes)):
            uniprot_acc = self.protein_nodes.loc[i, "uniprot_accession"]
            name = self.protein_nodes.loc[i, "uniprot_id"]
            self.names_mapping[uniprot_acc] = name

    def display_network(self):
        """
        Draw our graphical network 
        --------
        Params:
            - None: use init defined protein/drug values
        Returns:
            - Graphical Image representation
        """
        #create graph G based on edgelist
        G = nx.read_edgelist(self.network_edgelist)
        #get color indication for each node
        color = []
        for prot in list(G):
            indic = self.protein_nodes[self.protein_nodes['uniprot_accession']== prot]['indications'].values
            indic = indic[0]
            color.append(self.color_mapping[indic])
        plt.figure(figsize=(8,8), dpi=150) 
        #relabel from uniprot_accession to uniprot_id
        G = nx.relabel_nodes(G, self.names_mapping)
        #draw graph
        nx.draw(G, nodelist = list(G), with_labels = True, node_color = color)
        plt.savefig(self.output_file)



def main():
    # input variables
    network_edgelist = sys.argv[1]
    protein_nodes = sys.argv[2]
    output_file = sys.argv[3]

    network_ex = NetworkGraph(network_edgelist, protein_nodes, output_file)
    network_ex.display_network()

if __name__== "__main__":
    main()
