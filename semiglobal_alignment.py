import numpy as np
from Bio import pairwise2
from Bio import SeqIO
import sys
from numba import njit
from timer import Timer

directions = {0: (-1, 0),
              1: (-1, -1),
              2: (0, -1)}


def read_data(path="/home/krzysztof/Pobrane/rosalind_smgb.txt"):
    records = []
    for record in SeqIO.parse(path, "fasta"):
        records.append(record)
    return records


@njit
def fill_matrix(matrix, h_values, seq1, seq2):
    for i in range(1, matrix.shape[0]):
        for j in range(1, matrix.shape[1]):
            h_values[0][i][j] = matrix[i - 1][j] - 1
            h_values[2][i][j] = matrix[i][j - 1] - 1
            if seq1[i - 1] == seq2[j - 1]:
                h_values[1][i][j] = matrix[i - 1][j - 1] + 1
            else:
                h_values[1][i][j] = matrix[i - 1][j - 1] - 1
            matrix[i][j] = np.max(h_values[:, i, j])


def set_start_matrix(matrix, gap_penalty):
    matrix[0] = gap_penalty * np.arange(0, matrix.shape[1])
    matrix[:, 0] = gap_penalty * np.arange(0, matrix.shape[0])


def set_start_h_values(h_values):
    h_values[2][0] = np.ones((h_values.shape[2]))


def get_one_of_the_best(matrix, h_values, seq1, seq2):
    alignment1 = ""
    alignment2 = ""
    index_i, index_j = matrix.shape
    index_i -= 1
    index_j -= 1
    index_seq1, index_seq2 = -1, -1
    path = np.argmax(h_values, axis=0)
    while not index_i == index_j == 0:
        choice = path[index_i][index_j]
        if choice == 1:
            alignment1 = seq1[index_seq1] + alignment1
            alignment2 = seq2[index_seq2] + alignment2
            index_seq1 += directions[choice][0]
            index_seq2 += directions[choice][1]
            index_i += directions[choice][0]
            index_j += directions[choice][1]
        elif choice == 0:
            alignment1 = seq1[index_seq1] + alignment1
            alignment2 = '-' + alignment2
            index_seq1 += directions[choice][0]
            index_i += directions[choice][0]

        else:
            alignment1 = '-' + alignment1
            alignment2 = seq2[index_seq2] + alignment2
            index_j += directions[choice][1]
            index_seq2 += directions[choice][1]
    return alignment1, alignment2


try:
    DNA_s, DNA_t = read_data()
except ValueError:
    print("\nBłąd odczytu danych!\n"
          "Możliwość znalezienia dopasowania tylko między dwoma sekwencjami! "
          "Podaj prawidłowy format danych.")
    sys.exit(0)


seq_t = "CAGCACTTGGATTCTCGG" #str(DNA_t.seq)
seq_s = "CAGCGTGG" #str(DNA_s.seq)
gap_penalty = -1
matrix = np.zeros((len(seq_t) + 1, len(seq_s) + 1), dtype=np.int32)
h_values = np.zeros((3, len(seq_t) + 1, len(seq_s) + 1), dtype=np.int32)

with Timer() as t:
    set_start_matrix(matrix, gap_penalty)
    set_start_h_values(h_values)
    fill_matrix(matrix, h_values, seq_t, seq_s)
    alignment = get_one_of_the_best(matrix, h_values, seq_t, seq_s)


print('Czas zadania wyniósł: %.03f sec.' % t.interval)
print(alignment[0])
print(alignment[1])

#print(pairwise2.align.globalxx(human_str, rat_str))