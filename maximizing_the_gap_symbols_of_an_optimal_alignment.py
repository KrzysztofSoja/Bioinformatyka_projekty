from Bio import SeqIO
import numpy as np
from numba import njit
from timer import Timer


def read_data(path="/home/krzysztof/Pobrane/rosalind_mgap.txt"):
    records = []
    for record in SeqIO.parse(path, "fasta"):
        records.append(record)
    return records


@njit
def fill_matrix(matrix, seq1, seq2):
    h_values = np.zeros((3, matrix.shape[0], matrix.shape[1]))
    for i in range(1, matrix.shape[0]):
        for j in range(1, matrix.shape[1]):
            h_values[0][i][j] = matrix[i - 1][j] - 1
            h_values[1][i][j] = matrix[i][j - 1] - 1
            if seq1[i - 1] == seq2[j - 1]:
                h_values[2][i][j] = matrix[i - 1][j - 1] + 1
            else:
                h_values[2][i][j] = matrix[i - 1][j - 1] - 1
            matrix[i][j] = np.max(h_values[:, i, j])


def calculate_max_gap(matrix):
    return int(sum(matrix.shape) - 2*matrix[-1][-1])


with Timer() as t:
    seq1, seq2 = read_data()
    seq1, seq2 = str(seq1.seq), str(seq2.seq)
    seq1, seq2 = "AACGTA", "ACACCTA"
    alignment_matrix = np.zeros((len(seq1), len(seq2)))
    alignment_matrix[0] = -np.arange(0, alignment_matrix.shape[1])
    alignment_matrix[:, 0] = -np.arange(0, alignment_matrix.shape[0])
    fill_matrix(alignment_matrix, seq1, seq2)
    max_gap = calculate_max_gap(alignment_matrix)
print("Zadanie wykonane w czasie: %.3f sekund" % t.interval)
print(alignment_matrix)
print(max_gap)
