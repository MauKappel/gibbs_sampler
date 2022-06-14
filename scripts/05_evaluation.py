import numpy as np
# import matplotlib.pyplot as plt
from time import time
from sklearn.metrics import roc_curve, auc, matthews_corrcoef
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
from argparse import ArgumentParser
import pandas as pd


def plot_roc_curve():
    plt.title('Receiver Operating Characteristic')
    plt.plot(fpr, tpr, label='AUC = %0.3f (%smer)' % (roc_auc))
    plt.legend(loc='lower right')
    plt.plot([0, 1], [0, 1], c='black', linestyle='--')
    plt.ylabel('True Positive Rate')
    plt.xlabel('False Positive Rate')


eval_targets = np.loadtxt(eval_targets_file, dtype=str)
eval_peptides = np.loadtxt(eval_peptidesfile, dtype=str)

eval_targets_class = np.where(eval_targets > 0.426, 1, 0)
eval_peptides_class = np.where(eval_peptides > 0.426, 1, 0)
# Combining targets and prediction values with peptide length in a dataframe

plt.figure(figsize=(7, 7))
# For each peptide length compute AUC and plot ROC
fpr, tpr, threshold = roc_curve(eval_targets_class, eval_peptides_class)
roc_auc = auc(fpr, tpr)
plot_roc_curve()

mcc = matthews_corrcoef(eval_targets_class, eval_peptides_class)


def plot_mcc():
    plt.title('Matthews Correlation Coefficient')
    plt.scatter(eval_targets, eval_peptides, label='MCC = %0.3f' % mcc)
    plt.legend(loc='lower right')
    plt.ylabel('Predicted')
    plt.xlabel('Validation targets')
    plt.show()


plot_mcc()

pcc = pearsonr(eval_targets, eval_peptides)
print("PCC: ", pcc[0])

plt.scatter(eval_peptides, eval_targets);
