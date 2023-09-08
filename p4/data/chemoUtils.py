import numpy as np
import pandas as pd

def load_data(drug_file, target_file): 
    """
    Reads parameters from input files and stores them within the instance

    Inputs:
       drug_file = file with name, id, and fingerprint
       target_file = drug file with its protein target
    Outputs:
        self.drugs = data frame containing drug name/id and fingerprints
        self.targets = data frame containing drugs name/id and their protein targets
        
    """
    drugs = pd.read_csv(drug_file, sep=',')
    targets = pd.read_csv(target_file, sep = ',')
    return drugs, targets

def get_sets(drugs, targets, calc_targets):
    """
    Create drug and target set dictionaries for tanimoto analysis

    Inputs:
       drugs = dictionary for drugs and their fingerprints
       targets = dictionary for drugs and their protein targets
    Outputs:
        self.drug_setmap = dictionary with drugs as the keys and fingerprint sets as the values
        self.targ_setmap = dictionary with drugs as the keys and a set of protein targets as the values.
        
    """
    drug_setmap = {}
    targ_setmap = {}
    for i in np.arange(len(drugs)):
        drug = drugs['db_id'].loc[i]
        #parse fingerprint from dataframe and store in a set for efficient computation
        fingerprint = set([int(p) for p in drugs['maccs'].loc[i].split(' ')])
        drug_setmap[drug] = fingerprint
        #if necessary, store proteins in a set 
        if calc_targets:
            #get all proteins that bind a particular ligand and store in a set
            proteins = set(targets[targets['db_id'] == drug]['uniprot_accession'].values)
            targ_setmap[drug] = proteins
    return drug_setmap, targ_setmap
    

def get_tanimoto_val(fingerprint_drug_i, fingerprint_drug_j):
    """
    Get the tanimoto value based on fingerprint sets of both drugs.
    Tc = |fpt(mol_A) & fpt(mol_B)|/|fpt(mol_A) U fpt(mol_B)|

    Inputs:
       fingerprint_drug_i = fingerprint of drug i
       fingerprint_drug_j = fingerprint of drug j
    Outputs:
        tanimoto = tanimoto value calculated based on the Tc formula above
        
    """
    #intersection of drug fingerprints
    num = fingerprint_drug_i & fingerprint_drug_j 
    #union of drug fingerprints
    denom = fingerprint_drug_i | fingerprint_drug_j
    tanimoto = np.round(len(num)/len(denom), 6)
    return tanimoto