#!/usr/bin/env python3
# coding: utf-8

## IMPORT LIBRARIES

import os
import random
import subprocess
from joblib import Parallel, delayed
from time import time
from argparse import ArgumentParser

## PARSER
parser = ArgumentParser(description="Cross-validation")
parser.add_argument("-s", action="store", dest="seed", type=int, default=1, help="Random number seed (default: 1)")
parser.add_argument("-a", action="store", dest="allele", type=str, help="Allele (shared string)")
parser.add_argument("-k", action="store", dest="k", type=int, default=5, help="Cross-validation fold (default: 5)")
parser.add_argument("-m", action="store", dest="method", type=str, default="simple_gibbs", help="Method used in the algorithm (options: simple_gibbs(default), hobohm2, seq_weight)")

args = parser.parse_args()
seed = args.seed
allele = args.allele
k = args.k
method = args.method

## FUNCTIONS
def unix_call(command):
    '''Run jobs (unix command)'''
    job = subprocess.run(command.split())

## MAIN CODE

# Data path
script_path = os.getcwd()
data_dir = os.path.join(script_path, "../data/non-binders/")
results_dir = os.path.join(script_path, "../results/" + method)

# Define random seed
random.seed(seed)

# Retrieve evaluation file and run k CV fold Gibbs Sampler (1 method)
job = subprocess.run(["ls " + data_dir], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
infile_list = job.stdout.split("\n")
infile_list.pop()
if not os.path.exists(results_dir):
    os.makedirs(results_dir)
job_list = []
allele_list = []
for full_allele in infile_list:
    if allele in full_allele:
        allele_list.append(full_allele)
        if not os.path.exists(results_dir + "/" + full_allele):
            os.makedirs(results_dir + "/" + full_allele)
        if method == "hobohm2":
            hobohm_dir = os.path.join(script_path, "../data/hobohm2/" + full_allele)
            #os.makedirs(hobohm_dir)
            job = subprocess.run(["python3 02_hobohm2.py -f " + data_dir + full_allele + "/all" + " -d " + hobohm_dir +
                                  " -debug True -t 0.2"], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
            allele_dir = hobohm_dir + "/c00"
            evaluation_file = allele_dir + "0"
        else:
            allele_dir = data_dir + full_allele + "/c00"
            evaluation_file = allele_dir + "0"
        for i in range(1, k):
            training_index = list(range(1,k))
            training_index.remove(i)
            test_file = allele_dir + str(i)
            training_file = results_dir + "/" + full_allele + "/training_file_" + str(i)
            out_file_kld = results_dir + "/" + full_allele + "/kld_file_" + str(i)
            out_file_mat = results_dir + "/" + full_allele + "/mat_file_" + str(i)
            out_file_test = results_dir + "/" + full_allele + "/test_out_" + str(i)
            out_file_eval = results_dir + "/" + full_allele + "/eval_out_" + str(i)
            subprocess.run(["cat " + allele_dir + str(training_index[0]) + " " + allele_dir + str(training_index[1]) +
                           " " + allele_dir + str(training_index[2]) + " > " + training_file], shell=True,
                           stdout=subprocess.PIPE, universal_newlines=True)
            if method == "seq_weight":
                job_list.append("python3 03_gibbs_sampler.py -w -f1 " + training_file + " -f2 " + test_file + " -f3 " + evaluation_file + " -o1 " + out_file_kld + " -o2 " + out_file_mat + " -o3 " + out_file_test + " -o4 " + out_file_eval)
            else:
                job_list.append("python3 03_gibbs_sampler.py -f1 " + training_file + " -f2 " + test_file + " -f3 " + evaluation_file + " -o1 " + out_file_kld + " -o2 " + out_file_mat + " -o3 " + out_file_test + " -o4 " + out_file_eval)

# Parallelize Gibbs Sampler call
result = Parallel(n_jobs=8)(delayed(unix_call)(job) for job in job_list)

# Retrieve mean of all 4 evaluations
eval_job_list = []
for allele in allele_list:
    allele_dir = results_dir + "/" + allele + "/"
    eval_file = allele_dir + "eval_out_"
    mat_file = allele_dir + "mat_file_"
    output_file = allele_dir + "final_out"
    output_plot = allele_dir + "plot"
    eval_job_list.append("python3 05_evaluation_v2.py -k " + k + " -e " +  eval_file + " -mat " + mat_file + " -of " + output_file + " -op " + output_plot)
