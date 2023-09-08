
import sys
import numpy as np
import pandas as pd
import math


class GSEA(object):
    """
    GSEA class runs the GSEA algorithm by calculating enrichment scores and getting significant sets
    """
    def __init__(self):
        """
           Function to initialize objects holding gene info necessary for GSEA.
        """
        self.expression = pd.DataFrame()
        self.sample = pd.DataFrame()
        self.genesets = {}
        self.genes = []
        self.sample_lst = []
        self.all_dist = {}
        self.permuted_sample = False

    def load_data(self, expfile, sampfile, genesets): 
        """
        Reads parameters from input files and stores them within the GSEA instance

        Inputs:
           expression_file = file with expression levels for thousands of genes per sample
           sample_file = file listing each sample and whether it has the condition of interest
           genesets = file with geneset ids as the index and genes contained in the set as the value
        Outputs:
            self.expression = data frame containing expression level information per sample
            self.sample = sample data with corresponding boolean case/control column
            self.genesets = geneset data with set_id as key and gene list as value
            self.sample_lst = list of all sample ids
            self.genes = list of genes with calculated expression levels
        """
        self.expression = pd.read_csv(expfile, sep='\t', index_col=0)
        self.genes = self.expression.index
        self.sample = pd.read_csv(sampfile, sep='\t', header = None).rename(columns={1: 'score', 0: 'sample'}).set_index('sample')
        self.sample_lst =  self.sample.index
        var = True
        with open(genesets, 'r') as sets:
            while var:
                row = [i for i in sets.readline().strip().split('\t')]
                if row != ['']:
                    self.genesets[row[0]] = row[2:]
                else:
                    var = False



    def get_permutation(self):
        """
        Permutes the sample dataframe
        Outputs:
           sample_copy = permuted sample dataframe where the case/control booleans have been randomly 
           assigned to samples.
        """
        sample_copy = self.sample.copy()
        sample_copy['score'] = np.random.permutation(sample_copy['score'])
        return sample_copy

    
    def get_permuted_sample(self):
        """
        Returns the proper sample dataframe

        Outputs:
           Return a permuted sample if permuted_sample bool == True. 
           Otherwise, return the original sample dataframe,
        """
        if self.permuted_sample:
            return self.get_permutation()
        else:
            return self.sample


    def get_gene_rank_order(self):
        """
        Returns a list of all genes ranked by their logFC between patient and control

        Outputs:
           Returns a list of all genes (as strings) ranked by their logFC between patient and control, 
           with the gene with the highest logFC ordered first.
        """
        sample = self.get_permuted_sample()
        healthy = sample[sample['score'] == 0].index.tolist()
        patient = sample[sample['score'] == 1].index.tolist()
        healthy_gene_mean = np.mean(self.expression.loc[:, healthy], axis = 1)
        patient_gene_mean = np.mean(self.expression.loc[:, patient], axis = 1)
        logFCgene = (patient_gene_mean - healthy_gene_mean).tolist()
        logdiff_df = pd.DataFrame(data = {'gene': self.genes, 'logFC': logFCgene}).set_index(['gene'])
        logdiff_df = logdiff_df.sort_values('logFC', ascending = False)
        return logdiff_df.index.tolist()
    
    def get_enrichment_score(self, geneset):
        """
        Returns the enrichment score for a geneset using the GSEA random walk algorithm

        Inputs:
           geneset = string name of a geneset  
        Outputs:
           enrichment score (a float correct to two decimal places) for a given gene set
        """
        curr_genes = list(set(self.genesets[geneset]) & set(self.genes))
        G = len(curr_genes)
        N = len(self.genes)
        rank = self.get_gene_rank_order()
        running_sum = 0
        brownian_max = 0
        for i in rank:
            if i in curr_genes:
                running_sum += np.sqrt((N-G)/G)
            else:
                running_sum -= np.sqrt(G/(N-G))
            if running_sum > brownian_max:
                brownian_max = running_sum
        return np.round(float(brownian_max), 2)

    def get_sig_sets(self, p):
        """
        Returns the significant gene sets at a corrected threshold of p.

        Inputs:
           p = threshold parameter for determining significance
        Outputs:
           significant_sets = a list of significant gene sets (string, by name) at a corrected threshold of p. 
           If no gene sets are significant, significant_sets is an empty list.
        """
        significant_sets = []
        bonferroni_correction = p/len(self.genesets)
        for geneset in self.genesets:
            larger_count = 0
            actual_enrichment_score = self.get_enrichment_score(geneset)
            self.permuted_sample = True
            for i in np.arange(100):
                enrichment_score = self.get_enrichment_score(geneset)
                if enrichment_score >= actual_enrichment_score:
                    larger_count += 1
            if larger_count/100 <= bonferroni_correction:
                significant_sets.append(geneset)
            self.permuted_sample = False
        return significant_sets


def main():

    # check that the file is being properly used
    if (len(sys.argv) != 4):
        print("Please specify an input file and an output file as args.")
        return  
    # input variables
    expression_file = sys.argv[1]
    sample_file = sys.argv[2]
    KEGG_file = sys.argv[3]

    #testing a GSEA object and running 
    GSEA_ex = GSEA()
    GSEA_ex.load_data(expression_file, sample_file, KEGG_file)
    GSEA_ex.get_gene_rank_order()
    GSEA_ex.get_sig_sets(0.1)

if __name__=="__main__":
    main()
