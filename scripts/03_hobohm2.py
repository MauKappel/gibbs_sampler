import numpy as np
from time import time
import os
import argparse

parser = argparse.ArgumentParser(description="Hobohm2")
parser.add_argument("-f", action="store", dest="ALIGNEMENT_FILE", type=str, default="DRB1_0101/all", help="Alignement file, place in raw_data/Hobholm")
args = parser.parse_args()

script_path = os.getcwd()
data_dir = os.path.join(script_path, '../non-binders/')

# reading file
alignment_file_part = os.path.join(data_dir, f"{args.ALIGNEMENT_FILE}")
alignment_output = np.loadtxt(alignment_file_part, dtype=str)

print(alignment_output)