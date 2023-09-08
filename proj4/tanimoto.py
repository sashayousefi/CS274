import sys
import numpy as np
import pandas as pd
import chemoUtils
import itertools
import matplotlib.pyplot as plt

class Tanimoto(object):
    """
    Tanimoto class find the the Tanimoto Coefficient which is a metric to compare two chemicals 
    """
    def __init__(self, drug_file, target_file, output_file):
        """
            Function to initialize objects holding drug/protein info necessary for Tanimoto.
        """
        self.output_file = output_file
        self.drugs, self.targets = chemoUtils.load_data(drug_file, target_file) #load drug and target data into dataframes

        #fingerprint : dictionary of drug(key) and it's chemical fingerprint(value) (fingerprint = a set of 'keys' indicating the presence of a feature in the molecule)
        #protein : dictionary of drug(key) and proteins that bind to it (value)
        self.fingerprint, self.proteins = chemoUtils.get_sets(self.drugs, self.targets, True)

        self.tanimoto = {} #initialize a dictionary to hold all drug pairs and their tanimoto value
        self.targets_dict = {} #initialize a dictionary to indicate whether a pair of drugs shares a target
        self.output_df = pd.DataFrame()


    def calc_tanimoto(self):
        """
        Calculates the Tanimoto value (Jaccard Index)
        --------
        Params:
            - None: uses class fingerprint dictionaries
        Returns:
            - updates targets_dict: indicates whether a pair of drugs share a target
            - updates taninmoto_dict: indicates tanimoto value for a pair of drugs
        """
        for pair in itertools.combinations(self.drugs['db_id'], 2):
            drug_i = pair[0]
            drug_j = pair[1]
            #get the tanimoto values fo the drug pairs
            self.tanimoto[(drug_i, drug_j)] = chemoUtils.get_tanimoto_val(self.fingerprint[drug_i], \
            self.fingerprint[drug_j])
            targ1 = self.proteins[drug_i]
            targ2 = self.proteins[drug_j]
            #if drug sets share targets, set bool to 1, else set bool to 0 
            if (targ1 & targ2):
                self.targets_dict[pair] = 1
            else:
                self.targets_dict[pair] = 0
        
        
    def get_output(self):
        """
        Get the output dataframe of drug pairs, tanimoto values, and target boolean
        --------
        Params:
            - None: usees tanimoto dictionary
        Returns:
            - creates output dataframe with drug pairs, tanimoto values, and bool of whether the drugs share
              a target
        """
        #run calc_tanimoto function to get info for dataframe
        self.calc_tanimoto()
        first_drug = [i[0] for i in self.tanimoto.keys()]
        second_drug = [i[1] for i in self.tanimoto.keys()]
        #create output datqaframe
        self.output_df = pd.DataFrame(data = {'first': first_drug, 'second': second_drug, \
        'tanimoto':self.tanimoto.values(), 'shared': self.targets_dict.values()})
        self.output_df.to_csv(self.output_file, header=False, index = False)
    
    def create_images(self, data, name, file):
        """
        Create Tanimoto histogram using the Freedman-Diaconis bin width formula
        --------
        Params:
            - None
        Returns:
            - creates output histogram figure of tanimoto values. 
        """
        plt.clf() #clear plot

        #creating bins with freedman-diaconis bin width formula
        num = max(data) - min(data) #get data range
        quant3, quant1 = np.percentile(data, [75 ,25]) #get quantiles
        iqr = quant3 - quant1
        denom = 2*iqr/(len(data)**(1. / 3))

        #plot histogram of tanimoto values
        plt.hist(data, bins=int(num/denom))
        plt.title('syousefi {}'.format(name))
        plt.xlabel('Value bins')
        plt.ylabel('Number of occurences')
        plt.savefig(file)
    
    def run_images(self):
        """
        Plot Tanimoto histograms
        --------
        Params:
            - None: usees tanimoto dictionary and targets dictionary
        Returns:
            - all_tanimoto.png : histogram image of ALL tanimoto values
            - shared_tanimoto.png : histogram image of tanimoto values with SHARED targets
            - notshared_tanimoto.png : histogram image of tanimoto values with NO SHARED targets
        """
        pics =  ['all_tanimoto.png', 'shared_tanimoto.png', 'notshared_tanimoto.png']
        name = ['All', 'Shared', 'Not Shared']
        for i in np.arange(len(pics)):
            if i == 0:
                #all tanimoto vals
                data = list(self.tanimoto.values())
                self.create_images(data, name[i], pics[i])
            elif i == 1:
                #shared tanimoto vals
                pairs = [k for k,v in self.targets_dict.items() if v == 1]
                data = [self.tanimoto[k] for k in pairs]
                self.create_images(data, name[i], pics[i])
            else:
                #not shared tanimoto vals
                pairs = [k for k,v in self.targets_dict.items() if v == 0]
                data = [self.tanimoto[k] for k in pairs]
                self.create_images(data, name[i], pics[i])

def main():
    #input variables
    drug_file = sys.argv[1]
    target_file = sys.argv[2]
    output_file = sys.argv[3]

    #creating a Tanimoto object and run the algorithm/image generator    
    tanimoto_ex = Tanimoto(drug_file, target_file, output_file)
    tanimoto_ex.get_output()
    tanimoto_ex.run_images()


if __name__== "__main__":
    main()
