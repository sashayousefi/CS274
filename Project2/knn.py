import sys
import numpy as np
import pandas as pd


class KNN(object):
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

    def load_data(self, expfile, sampfile): 
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
        self.expression = pd.read_csv(expfile, sep='\t', index_col=0)
        self.sample = pd.read_csv(sampfile, sep='\t', index_col=0, header = None)
        self.sample =  pd.read_csv(sampfile, sep='\t', header = None).rename(columns={1: 'score', 0: 'sample'}).set_index('sample')
        self.sample_lst =  self.sample.index

    def calculate_distance(self, samp1, samp2):
        """
        Calculates the euclidian distance between two samples

        Inputs:
           samples = sample ids 
        Output:
            euclidian distance between the expression levels of sample1 and sample2. 
        """
        pt1 = self.expression.loc[:, samp1].to_numpy()
        pt2 = self.expression.loc[:, samp2].to_numpy()
        return np.linalg.norm(pt1 - pt2)


    def get_nearest_neighbors(self, samp1):     
        """
        Helper function that returns a sorted list of nearest neighbors for a given sample

        Inputs:
           sample = sample id 
        Output:
            nn = sorted list of sample ids where the nearest sample is first and the furthest 
            sample is last 
        """ 
        distances = []  
        for samp2 in self.sample_lst:
            if samp1 != samp2:
                if (samp2, samp1) not in self.all_dist.keys():
                    dist = self.calculate_distance(samp1, samp2)
                    self.all_dist[(samp1, samp2)] = dist
                else:
                    dist = self.all_dist[(samp2, samp1)]
                distances.append((samp2, dist))
        distances.sort(key=lambda tup: tup[1])
        nn = [dist[0] for dist in distances]
        return nn
        
    def get_assignments(self, k, fn):
        """
        Function that runs the KNN algorithm. Uses LOOCV to predict if each sample 
        is a case/control based on its k nearest neighbors.

        Inputs:
           k = number of nearest neighbors to consider in the algorithm
           fn = threshold parameter to distinguish between case & control 
        Output:
            list of the predicted case/control status for each sanmple. 
        """ 
        predicted_df = pd.DataFrame(columns = ["sample", "score"]).set_index("sample")
        nearest_neighbors = {}
        #LOOCV
        for test_case in self.sample_lst:
            nearest_neighbors[test_case] = self.get_nearest_neighbors(test_case)[:k]
            y_vals_train = [self.sample.loc[i] for i in nearest_neighbors[test_case]]
            if np.sum(y_vals_train) > k * fn:
                predicted_df.loc[test_case] = 1
            else:
                predicted_df.loc[test_case] = 0
        return list(predicted_df.loc[:, 'score'])

    def calc_metrics(self, k, fn):
        """
        Calculates accuracy metrics (sensitivity, specificity) for the KNN results

        Inputs:
           k = number of nearest neighbors to consider in the algorithm
           fn = threshold parameter to distinguish between case & control 
        Output:
            sensitivity = true positive rate 
            specificity = true negative rate
        """ 
        predictions = self.get_assignments(k, fn)
        zip_vals = list(zip(predictions, list(self.sample.loc[:, 'score'])))
        tp = np.sum([1 for i in zip_vals if i == (1, 1)])
        fp = np.sum([1 for i in zip_vals if i == (1, 0)])
        tn = np.sum([1 for i in zip_vals if i == (0, 0)])
        fn = np.sum([1 for i in zip_vals if i == (0, 1)])
        sensitivity = float(tp/(tp + fn))
        specificity = float(tn/(tn + fp))
        return [sensitivity, specificity]

def main():

    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    expression_file = sys.argv[1]
    sample_file = sys.argv[2]

    #testing a KNN object and running    
    KNN_ex = KNN()
    KNN_ex.load_data(expression_file, sample_file)
    KNN_ex.calc_metrics(6, 0.5)

if __name__== "__main__":
    main()
