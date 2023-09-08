import sys
import numpy as np
import pandas as pd
import itertools
import argparse
import chemoUtils


class pvalue(object):
    """
    P_values class to find the p_value for a observing a particular Tsummary value between two proteins.
    """
    def __init__(self, n, r, drug_file, target_file, proteinA, proteinB):
        """
            Function to initialize objects holding the necessary drug/protein info for p-value calculation
        """
        self.n = n #num iterations for bootstrapping
        self.r = r #random seed
        np.random.seed(self.r)
        self.proteinA = proteinA
        self.proteinB = proteinB
        self.drugs, self.targets = chemoUtils.load_data(drug_file, target_file)#load drug and target data into dataframes

        #set of ligands (drugs) that bind protein A and Protein B
        self.ligand_set_a = set(self.targets[self.targets['uniprot_accession'] == self.proteinA]['db_id'].values)
        self.ligand_set_b = set(self.targets[self.targets['uniprot_accession'] == self.proteinB]['db_id'].values)

        self.len_set_a, self.len_set_b = len(self.ligand_set_a), len(self.ligand_set_b)
        self.ligands = np.unique(self.drugs['db_id'].values) #unique ligand values to remove repeats
        #fingerprint : dictionary of drug(key) and it's chemical fingerprint(value) (fingerprint = a set of 'keys' indicating the presence of a feature in the molecule)
        self.fingerprints, _ = chemoUtils.get_sets(self.drugs, self.targets, False)
        self.tanimoto = {} #initialize a dictionary to hold all drug pairs and their tanimoto value

    
    def compute_T_summary(self, drug_pairs):
        """
        Compute the T_summary statistics : sum(Tc(a, b) | Tc > 0.5 for a in ligand_set_prot_A and b in ligand_set_prot_B)
        --------
        Params:
            - drug_pairs: all pairs of ligands between the set of ligands for protein A and the set of ligands for protein B
        Returns:
            - T_summary: the statistic calculated using the formula shown in this description
        """
        t_summary = 0
        for pairs in drug_pairs:
            drug_i = pairs[0]
            drug_j = pairs[1]
            #if the drugs are the same, simply return 1
            if drug_i == drug_j:
                tanimoto = 1
            else:
                #get cached results for tanimoto values
                if (drug_i, drug_j) in self.tanimoto.keys():
                    tanimoto = self.tanimoto[(drug_i, drug_j)]
                else:
                    #cache drug pair and taminoto value in a dictionary for faster computation
                    tanimoto = chemoUtils.get_tanimoto_val(self.fingerprints[drug_i], self.fingerprints[drug_j])
                    self.tanimoto[(drug_i, drug_j)] = tanimoto
            #only add to sum if tanimoto > 0.5
            if tanimoto > 0.5:
                t_summary += tanimoto
        return t_summary


    def bootstrap(self, print_pval = True):
        """
        Calculate p_value using bootstrap sampling.
        This bootstrap sampling procedure will give us an empirical estimate of 
        the probability of getting the value of Tsummary(A,B) by selecting two proteins at random.
        Params:
            - print_pval: bool - should we print the p value
        Returns:
            - p_bootstrap: the p-value attained from the bootstrap sampling procedure (how likely are we to 
              have observed this Tsummary(A, B) value with random chance)
        """
        Tb_summary = np.array([])
        #compute the T_summary(A,B) value by using the set of ligands for protein A and the set of ligands for protein B
        T_summary = self.compute_T_summary(itertools.product(self.ligand_set_a, self.ligand_set_b))
        #iterate for n number of bootstraps
        for i in np.arange(self.n):
            #get random sets of ligands for protein A and B which will be used to determine if our T_summary(A, B)
            #was likely to come from random chance.
            random_sets_A = np.random.choice(self.ligands, size=self.len_set_a, replace=True)
            random_sets_B = np.random.choice(self.ligands, size=self.len_set_b, replace=True)
            #get all pairs of ligands between the two random sets
            drug_pairs = itertools.product(random_sets_A, random_sets_B)   
            #comput T_summary for the random sets
            Tb_summary = np.append(Tb_summary, self.compute_T_summary(drug_pairs))
        #calculate p_value using T_summary values from bootstrapping procedure
        p_bootstrap = np.sum(Tb_summary > T_summary)/self.n
        if print_pval:
            print(p_bootstrap)
        return p_bootstrap


def main():   
    def get_args():
        """	
	    Function to get the arguments from the command line input. Takes in a iteration number (default 500), random seed (default 214),
	    drug file, target file, protein A, and protein B.

        Returns: 
           - Init arguments
	    """
        parser = argparse.ArgumentParser(description='Read in optional info')
        parser.add_argument('-n', help='the number of iterations', metavar=int, default = 500)
        parser.add_argument('-r', help='parameter that sets the state of the pseudo-random number generator in Python.', metavar=int, default = 214)
        parser.add_argument('drug_file')    
        parser.add_argument('target_file')    
        parser.add_argument('proteinA')   
        parser.add_argument('proteinB')
        return parser.parse_args()

    #input variables
    args = get_args()
    n, r, drug_file, target_file, proteinA, proteinB = int(args.n), int(args.r), args.drug_file, args.target_file, args.proteinA, args.proteinB
    pvalue_ex = pvalue(n, r, drug_file, target_file, proteinA, proteinB)
    pvalue_ex.bootstrap()

if __name__== "__main__":
    main()
