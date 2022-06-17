import numpy as np
# import matplotlib.pyplot as plt
from time import time
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import pandas as pd
import logomaker as lm
import os


parser = ArgumentParser(description="GibbsSampler evaluation")

parser.add_argument("-k", action="store", dest="k", type=int, default=4, help="Number of cross-validation folds (default: 4)")
parser.add_argument("-p", action="store", dest="path", type=str, help="Path with evaluation data and weighted scores")
parser.add_argument("-of", action="store", dest="output_file", type=str, help="Output file")
parser.add_argument("-op", action="store", dest="output_plot", type=str, help="Output plot file")

args = parser.parse_args()
k = args.k
path = args.path
output_file = args.output_file
output_plot = args.output_plot

def plot_mcc(k):
    plt.title('Matthews Correlation Coefficient')
    plt.scatter(eval_targets, eval_prediction, label='MCC = %0.3f' % mcc)
    plt.legend(loc='lower right')
    plt.ylabel('Predicted')
    plt.xlabel('Validation targets')
    plt.savefig(output_plot + "_mcc_" + str(k) + ".pdf")

def plot_roc_curve(k):
    plt.title('Receiver Operating Characteristic')
    plt.plot(fpr, tpr, label='AUC = %0.3f' % (roc_auc))
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], c='black', linestyle='--')
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig(output_plot + "_roc_" + str(k) + ".pdf")

def plot_pcc(k):
    plt.title('Pearson Correlation Coefficient')
    plt.scatter(eval_targets, eval_prediction, label='PCC = %0.3f' % pcc[0])
    plt.legend(loc='lower right')
    plt.ylabel('Predicted')
    plt.xlabel('Validation targets')
    plt.savefig(output_plot + "_pcc_" + str(k) + ".pdf")

def SeqPlot(weight_matrix):
    """
    Creates a simple sequence logo from a weight matrix
    Weight matrix is produced via 'from_psi_blast' function
    """
    SeqPlot = lm.Logo(df = pd.DataFrame(weight_matrix),
                    fade_below=0.5,
                    shade_below=0.5,
                    figsize=(15,7))

    # set axis labels
    SeqPlot.ax.set_xlabel('Amino acid position', fontsize=14)
    SeqPlot.ax.set_ylabel('bits', fontsize=14)

    # Highlighting specific or range of amino acids
    #SeqPlot.highlight_position(2,color='lightgray',alpha=0.5)
    #SeqPlot.highlight_position_range(2,6,alpha=0.5,color='lightgray')
    return SeqPlot

def initialize_matrix(peptide_length, alphabet):
    init_matrix = [0]*peptide_length

    for i in range(0, peptide_length):

        row = {}

        for letter in alphabet:
            row[letter] = 0.0

        init_matrix[i] = row

    return init_matrix

def from_psi_blast(file_name):
    f = open(file_name, "r")

    nline = 0
    for line in f:

        sline = str.split( line )

        if nline == 0:
        # recover alphabet
            alphabet = [str]*len(sline)
            for i in range(0, len(sline)):
                alphabet[i] = sline[i]

            matrix = initialize_matrix(peptide_length, alphabet)

        else:
            i = int(sline[0])

            for j in range(2,len(sline)):
                matrix[i-1][alphabet[j-2]] = float(sline[j])

        nline+= 1

    return matrix


for i in range(1, k):
    eval_data = np.loadtxt(path + "/eval_out_" + str(i), dtype=str).reshape(-1,5)

    eval_peptides = eval_data[1:, 2]
    eval_prediction = eval_data[1:, 3].astype(float)
    eval_targets = eval_data[1:, 4].astype(float)

    eval_targets_class = np.where(eval_targets > 0.426, 1, 0)
    eval_prediction_class = np.where(eval_prediction > 0.426, 1, 0)
    # Combining targets and prediction values with peptide length in a dataframe

    #------------------------------------------------------------------------------------#
    # ROC
    #------------------------------------------------------------------------------------#
    plt.figure(figsize=(7, 7))
    # For each peptide length compute AUC and plot ROC
    fpr, tpr, threshold = roc_curve(eval_targets_class, eval_prediction)
    roc_auc = auc(fpr, tpr)
    plot_roc_curve(i)

    #------------------------------------------------------------------------------------#
    # MCC
    #------------------------------------------------------------------------------------#
    plt.figure(figsize=(7, 7))
    mcc = matthews_corrcoef(eval_targets_class, eval_prediction_class)
    plot_mcc(i)

    #------------------------------------------------------------------------------------#
    # PCC
    #------------------------------------------------------------------------------------#
    plt.figure(figsize=(7, 7))
    pcc = pearsonr(eval_targets, eval_prediction)
    plot_pcc(i)

    #------------------------------------------------------------------------------------#
    # Sequence logo
    #------------------------------------------------------------------------------------#
    peptide_length = len(eval_peptides[0])
    #peptide_length = 9

    # Conversion from psi-blast to dictionary format
    w_matrix = from_psi_blast(path + "/mat_file_" + str(i))
    print(w_matrix, file = open(path + "/output_" + str(i), 'w'))

    SeqPlot(w_matrix)
    logo_name = path + "/seqlogo_" + str(i) + ".pdf"
    plt.savefig(logo_name)