import numpy as np
# import matplotlib.pyplot as plt
from time import time
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from argparse import ArgumentParser
from os.path import isdir, isfile
from os import listdir
import pandas as pd
import logomaker as lm
import os

parser = ArgumentParser(description="GibbsSampler evaluation")

parser.add_argument("-k", action="store", dest="k", type=int, default=5,
                    help="Number of cross-validation folds (default: 5)")
parser.add_argument("-p", action="store", dest="path", type=str,
                    help="Path with evaluation data/ weighted scores/ output path (e.g. simple_gibbs)")

args = parser.parse_args()
path = args.path
k = args.k


def plot_pcc(path_out):
    plt.title('Pearson Correlation Coefficient')
    plt.scatter(eval_targets, eval_prediction, label='PCC = %0.3f' % pcc[0])
    plt.legend(loc='lower right')
    plt.ylabel('Predicted')
    plt.xlabel('Validation targets')
    plt.savefig(path_out + "/pcc.pdf")


def plot_roc_curve(path_out):
    plt.title('Receiver Operating Characteristic')
    plt.plot(fpr, tpr, label='AUC = %0.3f' % roc_auc)
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], c='black', linestyle='--')
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')
    plt.savefig(path_out + "/roc.pdf")


def plot_mcc(path_out):
    plt.title('Matthews Correlation Coefficient')
    plt.scatter(eval_targets, eval_prediction, label='MCC = %0.3f' % mcc)
    plt.legend(loc='lower right')
    plt.ylabel('Predicted')
    plt.xlabel('Validation targets')
    plt.savefig(path_out + "/mcc.pdf")


eval_peptides = np.array([])
eval_prediction = []
eval_targets = []
eval_targets_class = []
eval_prediction_class = []

# Iterate through folders in path (e.g. results/simple_gibbs/DRB1_...)
for eval_file in list(filter(lambda x: isdir(f"{path}\\{x}"), listdir(path))):
    for i in range(1, k):
        eval_data = np.loadtxt(eval_file + "/eval_out_" + i, dtype=str).reshape(-1, 5)
        eval_peptides = np.concatenate((eval_peptides, eval_data[1:, 2]), axis=0)
        eval_prediction = eval_prediction.append(eval_data[1:, 3].astype(float))
        eval_targets = eval_targets.append(eval_data[1:, 4].astype(float))

        eval_targets_class = eval_targets.append(np.where(eval_targets > 0.426, 1, 0))
        eval_prediction_class = eval_prediction_class.append(np.where(eval_prediction > 0.426, 1, 0))
        # Combining targets and prediction values with peptide length in a dataframe

# ------------------------------------------------------------------------------------#
# ROC
# ------------------------------------------------------------------------------------#
plt.figure(figsize=(7, 7))
# For each peptide length compute AUC and plot ROC
fpr, tpr, threshold = roc_curve(eval_targets_class, eval_prediction)
roc_auc = auc(fpr, tpr)
plot_roc_curve(path)

# ------------------------------------------------------------------------------------#
# MCC
# ------------------------------------------------------------------------------------#
plt.figure(figsize=(7, 7))
mcc = matthews_corrcoef(eval_targets_class, eval_prediction_class)
plot_mcc(path)

# ------------------------------------------------------------------------------------#
# PCC
# ------------------------------------------------------------------------------------#
plt.figure(figsize=(7, 7))
pcc = pearsonr(eval_targets, eval_prediction)
plot_pcc(path)

print("Unique peptides:" + np.unique(eval_peptides).shape)
