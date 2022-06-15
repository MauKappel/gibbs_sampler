import numpy as np
import os
from argparse import ArgumentParser
import sys
import matplotlib.pylab as plt

parser = ArgumentParser(description="Alignment_O2, if only q or db is provided, sequences will be align within file")
parser.add_argument("-q", action="store", dest="query_file", default=None, help="File with query sequence")
parser.add_argument("-db", action="store", dest="db_file", default=None,help="File with database sequence")
parser.add_argument("-go", action="store", dest="gap_open", type=float, default=-11.0, help="Value of gap open (-11.0)")
parser.add_argument("-debug", action="store", dest="DEBUG", type=bool, default=False, help="only align a small number of sequences")
parser.add_argument("-ge", action="store", dest="gap_extension", type=float, default=-1.0,
                    help="Value of gap extension (-1.0)")
parser.add_argument("-t", action="store", dest="thresh", type=float, default=0.8, help="Value of similarity theshold (0.8)")

args = parser.parse_args()
gap_open = args.gap_open
gap_extension = args.gap_extension

## If only one file is provided allign all vs all within the file, else align 2 files
if args.db_file == None and args.query_file == None:
    print('provide alignment file')
    sys.exit()

elif args.query_file == None:
    query_file = args.db_file
    database_file = args.db_file

elif args.db_file == None:
    database_file = args.query_file
    query_file = args.query_file

else:
    query_file = args.query_file
    database_file = args.db_file

script_path = os.getcwd()
data_dir = os.path.join(script_path, '../data/')

### Amino Acid codes
alphabet_file = data_dir + "Matrices/alphabet"
alphabet = np.loadtxt(alphabet_file, dtype=str)

### Blosum Matrix
blosum_file = data_dir + "Matrices/BLOSUM50"
_blosum50 = np.loadtxt(blosum_file, dtype=float).reshape((24, -1)).T

blosum50 = {}

for i, letter_1 in enumerate(alphabet):

    blosum50[letter_1] = {}

    for j, letter_2 in enumerate(alphabet):
        blosum50[letter_1][letter_2] = _blosum50[i, j]


## Alignment Matrix
# This functions returns, apart from the final Alignment Matrix, all the intermedite Matrices (for plotting purposes).
def smith_waterman_alignment(query="VLLP", database="VLILP", scoring_scheme={}, gap_open=-5, gap_extension=-1):
    # Matrix imensions
    M = len(query)
    N = len(database)

    # D matrix change to float
    D_matrix = np.zeros((M + 1, N + 1), np.int32)

    # P matrix
    P_matrix = np.zeros((M + 1, N + 1), np.int32)

    # Q matrix
    Q_matrix = np.zeros((M + 1, N + 1), np.int32)

    # E matrix
    E_matrix = np.zeros((M + 1, N + 1), dtype=object)

    # Initialize matrices
    for i in range(M, 0, -1):
        # alignment_matrix[i-1, N] = alignment_matrix[i, N] + gap_open
        D_matrix[i - 1, N] = 0
        P_matrix[i - 1, N] = 0
        Q_matrix[i - 1, N] = 0
        E_matrix[i - 1, N] = 0

    for j in range(N, 0, -1):
        # alignment_matrix[M, j-1] = alignment_matrix[M, j] + gap_open
        D_matrix[M, j - 1] = 0
        P_matrix[M, j - 1] = 0
        Q_matrix[M, j - 1] = 0
        E_matrix[M, j - 1] = 0

    # Main loop
    D_matrix_max_score, D_matrix_i_max, D_matrix_j_max = -9, -9, -9
    for i in range(M - 1, -1, -1):  # Why not to 0 only?
        for j in range(N - 1, -1, -1):

            # Q_matrix[i,j] entry
            # gap_open_database = XX
            # gap_extension_database = XX
            gap_open_database = D_matrix[i + 1, j] + gap_open
            gap_extension_database = Q_matrix[i + 1, j] + gap_extension
            max_gap_database = max(gap_open_database, gap_extension_database)

            Q_matrix[i, j] = max_gap_database

            # P_matrix[i,j] entry
            gap_open_query = D_matrix[i, j + 1] + gap_open
            # gap_extension_query = XX
            gap_extension_query = P_matrix[i, j + 1] + gap_extension
            max_gap_query = max(gap_open_query, gap_extension_query)

            P_matrix[i, j] = max_gap_query

            # D_matrix[i,j] entry
            # diagonal_score = XX
            diagonal_score = D_matrix[i + 1, j + 1] + scoring_scheme[query[i]][database[j]]

            # E_matrix[i,j] entry
            candidates = [(1, diagonal_score),
                          (2, gap_open_database),
                          (4, gap_open_query),
                          (3, gap_extension_database),
                          (5, gap_extension_query)]

            direction, max_score = max(candidates, key=lambda x: x[1])

            # check entry sign
            if max_score > 0:
                E_matrix[i, j] = direction
                D_matrix[i, j] = max_score
            else:
                E_matrix[i, j] = 0
                D_matrix[i, j] = 0

            # fetch global max score
            if max_score > D_matrix_max_score:
                D_matrix_max_score = max_score
                D_matrix_i_max = i
                D_matrix_j_max = j

    return P_matrix, Q_matrix, D_matrix, E_matrix, D_matrix_i_max, D_matrix_j_max, D_matrix_max_score


# ## Alignment Matrix Traceback
def smith_waterman_traceback(E_matrix, D_matrix, i_max, j_max, query="VLLP", database="VLILP", gap_open=-5,
                             gap_extension=-1):
    M = len(query)
    N = len(database)

    aligned_query = []
    aligned_database = []

    # start from max_i, max_j
    i, j = i_max, j_max
    matches = 0
    while i < M and j < N:

        # E[i,j] = 0, stop back tracking
        if E_matrix[i, j] == 0:
            break

        # E[i,j] = 1, match
        if E_matrix[i, j] == 1:
            aligned_query.append(query[i])
            aligned_database.append(database[j])
            if (query[i] == database[j]):
                matches += 1
            i += 1
            j += 1

        # E[i,j] = 2, gap opening in database
        if E_matrix[i, j] == 2:
            aligned_database.append("-")
            aligned_query.append(query[i])
            i += 1

        # E[i,j] = 3, gap extension in database
        if E_matrix[i, j] == 3:

            count = i + 2
            score = D_matrix[count, j] + gap_open + gap_extension

            # Find length of gap
            while ((score - D_matrix[i, j]) * (score - D_matrix[i, j]) >= 0.00001):
                count += 1
                # score = XX
                score = D_matrix[count, j] + gap_open + (count - i - 1) * gap_extension

            for k in range(i, count):
                aligned_database.append("-")
                aligned_query.append(query[i])
                i += 1

        # E[i,j] = 4, gap opening in query
        if E_matrix[i, j] == 4:
            aligned_query.append("-")
            aligned_database.append(database[j])
            j += 1

        # E[i,j] = 5, gap extension in query
        if E_matrix[i, j] == 5:

            count = j + 2
            score = D_matrix[i, count] + gap_open + gap_extension

            # Find length of gap
            while ((score - D_matrix[i, j]) * (score - D_matrix[i, j]) >= 0.0001):
                count += 1
                score = D_matrix[i, count] + gap_open + (count - j - 1) * gap_extension

            for k in range(j, count):
                aligned_query.append("-")
                aligned_database.append(database[j])
                j += 1

    return aligned_query, aligned_database, matches


### Align query against database#

if args.DEBUG == False:
    print(f'aligning {query_file} vs {database_file}')
    query_list = np.loadtxt(query_file, dtype=str, delimiter=' ').reshape(-1, 2)
    database_list = np.loadtxt(database_file, dtype=str).reshape(-1, 2)
else:
    print('debug_mode')
    query_list = np.loadtxt(query_file, dtype=str, delimiter=' ').reshape(-1, 2)[0:100]
    database_list = np.loadtxt(database_file, dtype=str).reshape(-1, 2)[0:100]


from time import time

scoring_scheme = blosum50

# this returns current timestamp in seconds
t0 = time()
all_matches = []
match_matrix = np.zeros((len(query_list), len(database_list)))
for m, query in enumerate(query_list):
    query_sequence = query[0]

    for n, database in enumerate(database_list):
        database_sequence = database[0]

        P_matrix, Q_matrix, D_matrix, E_matrix, i_max, j_max, max_score = smith_waterman_alignment(query_sequence,
                                                                                                   database_sequence,
                                                                                                   scoring_scheme,
                                                                                                   gap_open,
                                                                                                   gap_extension)
        aligned_query, aligned_database, matches = smith_waterman_traceback(E_matrix, D_matrix, i_max, j_max,
                                                                            query_sequence, database_sequence, gap_open,
                                                                            gap_extension)
        '''
        TODO:
        - change the code so it only accepts one file
        - calculate the sum of the rows, remove the rows and corresponding column which have the highest count
        
        NOTE: 
        - Theshold greater than 80% matches will be a 1, lower than 80% matches will be 0
        '''
        #print("ALN", query_sequence, len(query_sequence), database_sequence, len(database_sequence), len(aligned_query),
              #matches, max_score)
        #print("QAL", i_max, ''.join(aligned_query))
        #print("DAL", j_max, ''.join(aligned_database))
        if query_sequence != database_sequence:
            if matches > args.thresh * len(query_sequence):
                all_matches.append(matches)
                match_matrix[m][n] = 1
            else:
                all_matches.append(matches)
                match_matrix[m][n] = 0

t1 = time()
print("Time (m):,", (t1 - t0) / 60)
print(match_matrix)
print(np.apply_along_axis(sum, axis=1, arr=match_matrix))
plt.hist(all_matches)
plt.savefig('../results/02_alignement_matches.png')