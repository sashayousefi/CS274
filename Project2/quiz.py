import numpy as np
import sys
import pandas as pd
import matplotlib.pyplot as plt
from knn import KNN



class ROC(object):
    def __init__(self, ks, fn, expression_file, sample_file):
        self.dict = {}
        self.k = ks
        self.fns = fn
        self.expression_file = expression_file
        self.sample_file = sample_file

    def create_ROC_vals(self):
        knn_ex = KNN()
        knn_ex.load_data(self.expression_file, self.sample_file)
        for fn in self.fns:
            self.dict[fn] = knn_ex.calc_metrics(self.k, fn)

    def create_ROC_curve(self):
        self.create_ROC_vals()
        tpr = [i[0] for i in self.dict.values()]
        fpr = [(1 - i[1]) for i in self.dict.values()]
        plt.figure()
        lw = 2
        plt.plot(
            fpr,
            tpr,
            color="navy",
            lw=lw,
            label="ROC curve",
        )
        plt.xlim([0.0, 1.0])
        plt.ylim([0.0, 1.05])
        plt.xlabel("False Positive Rate")
        plt.ylabel("True Positive Rate")
        plt.title("ROC curve for KNN with K = 3 and fns = {}".format(self.fns))
        plt.legend(loc="lower right")
        plt.show()



def main():
    # check that the file is being properly used
    if (len(sys.argv) !=3):
        print("Please specify an input file and an output file as args.")
        return
        
    # input variables
    expression_file = sys.argv[1]
    sample_file = sys.argv[2]
    lst = [0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 1]
    ROC_ex = ROC(3, lst, expression_file, sample_file)
    ROC_ex.create_ROC_curve()

if __name__== "__main__":
    main()
