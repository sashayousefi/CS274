import sys
import numpy as np
import pandas as pd


class Tanimoto(object):
    """
    KNN class runs the K nearest neighbors algorithm and calculates accuracy metrics
    """
    def __init__(self):
        """
            Function to initialize objects holding gene info necessary for KNN.
        """
        self.expression = pd.DataFrame()
        self.sample = pd.DataFrame()
        self.sample_lst = []
        self.all_dist = {}

    def load_data(self, drug_file, target_file): 
        """
        Reads parameters from input files and stores them within the KNN instance

        Inputs:
           expression_file = file with expression levels for thousands of genes per sample
           sample_file = file listing each sample and whether it has the condition of interest
        Outputs:
            self.expression = data frame containing expression level information per sample
            self.sample = sample data with corresponding boolean case/control column
            self.sample_lst = list of all sample ids
        """
        drugs = pd.read_csv(drug_file, sep='\t')
        print(drugs)
        '''self.expression = pd.read_csv(expfile, sep='\t', index_col=0)
        self.sample = pd.read_csv(sampfile, sep='\t', index_col=0, header = None)
        self.sample =  pd.read_csv(sampfile, sep='\t', header = None).rename(columns={1: 'score', 0: 'sample'}).set_index('sample')
        self.sample_lst =  self.sample.index'''

def main():

    # check that the file is being properly used
    if (len(sys.argv) != 4):
        print("Incorrect number of inputs")
        return
        
    # input variables
    drug_file = sys.argv[1]
    target_file = sys.argv[2]
    outputfile = sys.argv[3]

    #testing a KNN object and running    
    tanimoto_ex = Tanimoto()
    tanimoto_ex.load_data(drug_file, target_file)
    #KNN_ex.calc_metrics(6, 0.5)

if __name__== "__main__":
    main()
