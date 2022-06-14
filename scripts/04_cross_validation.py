#!/usr/bin/env python3
# coding: utf-8

## IMPORT LIBRARIES

import os
import random
import numpy as np
import math
from argparse import ArgumentParser

## PARSER
parser = ArgumentParser(description="Cross-validation")
parser.add_argument("-s", action="store", dest="seed", type=int, default=1, help="Random number seed (default: 1)")

args = parser.parse_args()
seed = args.seed

## FUNCTIONS

## MAIN CODE

# Data path
script_path = os.getcwd()
data_dir = os.path.join(script_path, '../non-binders/')

# Define random seed
random.seed(seed)

# Retrieve evaluation file


# Define training files and test file
