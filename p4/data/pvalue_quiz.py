import sys
import numpy as np
import pandas as pd
import itertools
import argparse
import chemoUtils


class pvalue(object):
    """
    KNN class runs the K nearest neighbors algorithm and calculates accuracy metrics
    """
    def __init__(self, n, r, drug_file, target_file, proteinA, proteinB):
        """
            Function to initialize objects holding gene info necessary for KNN.
        """
        self.n = n
        self.r = r
        np.random.seed(self.r)
        self.proteinA = proteinA
        self.proteinB = proteinB

        self.drugs, self.targets = chemoUtils.load_data(drug_file, target_file)
        self.ligand_set_a = set(self.targets[self.targets['uniprot_accession'] == self.proteinA]['db_id'].values)
        self.ligand_set_b = set(self.targets[self.targets['uniprot_accession'] == self.proteinB]['db_id'].values)
        self.len_set_a, self.len_set_b = len(self.ligand_set_a), len(self.ligand_set_b)
        self.ligands = np.unique(self.drugs['db_id'].values)
        self.fingerprints, _ = chemoUtils.get_sets(self.drugs, self.targets, False)
        self.tanimoto = {}

    
    def compute_T_summary(self, drug_pairs):
        t_summary = 0
        for pairs in drug_pairs:
            drug_i = pairs[0]
            drug_j = pairs[1]
            if drug_i == drug_j:
                tanimoto = 1
            else:
                if (drug_i, drug_j) in self.tanimoto.keys():
                    tanimoto = self.tanimoto[(drug_i, drug_j)]
                else:
                    tanimoto = chemoUtils.get_tanimoto_val(self.fingerprints[drug_i], self.fingerprints[drug_j])
                    self.tanimoto[(drug_i, drug_j)] = tanimoto
            if tanimoto > 0.5:
                t_summary += tanimoto
        return t_summary


    def bootstrap(self, print_pval = True):
        Tb_summary = np.array([])
        T_summary = self.compute_T_summary(itertools.product(self.ligand_set_a, self.ligand_set_b))
        for i in np.arange(self.n):
            random_sets_A = np.random.choice(self.ligands, size=self.len_set_a, replace=True)
            random_sets_B = np.random.choice(self.ligands, size=self.len_set_b, replace=True)
            drug_pairs = itertools.product(random_sets_A, random_sets_B)   
            Tb_summary = np.append(Tb_summary, self.compute_T_summary(drug_pairs))
        p_bootstrap = np.sum(Tb_summary > T_summary)/self.n
        if print_pval:
            print(p_bootstrap)
        return p_bootstrap
        




def main():
        
    # input variables
    def get_args():
        parser = argparse.ArgumentParser(description='Read in optional info')
        parser.add_argument('-n', help='the number of iterations', metavar=int, default = 500)
        parser.add_argument('-r', help='parameter that sets the state of the pseudo-random number generator in Python.', metavar=int, default = 214)
        parser.add_argument('drug_file')    
        parser.add_argument('target_file')    
        parser.add_argument('proteinA')   
        parser.add_argument('proteinB')
        return parser.parse_args()

    '''args = get_args()
    n, r, drug_file, target_file, proteinA, proteinB = int(args.n), int(args.r), args.drug_file, args.target_file, args.proteinA, args.proteinB
    pvalue_ex = pvalue(n, r, drug_file, target_file, proteinA, proteinB)
    pvalue_ex.bootstrap()'''

        
    def run_quiz():
        f = open('output_quiz.txt', 'w')
        for iter in [100, 500, 1000]:
            vals = []
            for seed in np.arange(100):
                p_val = pvalue(iter, i, 'drugs.csv', 'targets.csv', 'P54577', 'Q7RTX0')
                val = p_val.bootstrap(print_pval = False)
                vals += [val]
                f.write(str(iter) + '\t' + str(val))
                f.write('\n')
            print(np.mean(vals))
            print(np.std(vals))
        
    
    run_quiz()


if __name__== "__main__":
    main()
