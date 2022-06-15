#!/usr/bin/env python3
# coding: utf-8

## IMPORT LIBRARIES

import os
import random
import subprocess
from argparse import ArgumentParser

## PARSER
parser = ArgumentParser(description="Cross-validation")
parser.add_argument("-s", action="store", dest="seed", type=int, default=1, help="Random number seed (default: 1)")
parser.add_argument("-a", action="store", dest="allele", type=str, help="Allele (shared string)")
parser.add_argument("-k", action="store", dest="k", type=int, default=5, help="Cross-validation fold (default: 5)")

args = parser.parse_args()
seed = args.seed
allele = args.allele
k = args.k

## FUNCTIONS

## MAIN CODE

# Data path
script_path = os.getcwd()
data_dir = os.path.join(script_path, '../data/non-binders/')

# Define random seed
random.seed(seed)

# Retrieve evaluation file and run k CV fold Gibbs Sampler
job = subprocess.run(['ls ' + data_dir], shell=True, stdout=subprocess.PIPE, universal_newlines=True)
infile_list = job.stdout.split('\n')
infile_list.pop()
for element in infile_list:
    if allele in element:
        allele_dir = data_dir + element + "/c00"
        evaluation_file = allele_dir + "0"
        for i in range(1, 5):
            training_index = [1, 2, 3, 4]
            training_index.remove(i)
            test_file = allele_dir + str(i)
            training_file = data_dir + element + "/training_file_" + str(i)
            out_file_name = data_dir + element + "/mat_file_" + str(i)
            subprocess.run(['cat ' + allele_dir + str(training_index[0]) + " " + allele_dir + str(training_index[1]) +
                            " " + allele_dir + str(training_index[2]) + " > " + training_file], shell=True,
                           stdout=subprocess.PIPE, universal_newlines=True)
            subprocess.run(['python3 03_gibbs_sampler.py -f ' + training_file + " -o " + out_file_name], shell=True,
                           stdout=subprocess.PIPE, universal_newlines=True)
