import sys
import numpy as np
from Bio import SeqIO
from numba import njit
from timer import Timer


def read_data(path="/home/krzysztof/Pobrane/test.txt"):
    records = []
    for record in SeqIO.parse(path, "fasta"):
        records.append(record)
    return records


@njit
def fill_matrix(matrix, seq1, seq2):
    h_values = np.zeros((3,))
    for i in range(1, matrix.shape[0]):
        for j in range(1, matrix.shape[1]):
            h_values[0] = matrix[i - 1][j] - 1
            h_values[1] = matrix[i][j - 1] - 1
            if seq1[i - 1] == seq2[j - 1]:
                h_values[2] = matrix[i - 1][j - 1] + 1
            else:
                h_values[2] = matrix[i - 1][j - 1] - 1
            matrix[i][j] = np.max(h_values)


#@njit
def danied_equal_fill_matrix(matrix, seq1, seq2, d_index_i, d_index_j):
    h_values = np.zeros((3,))
    for i in range(1, matrix.shape[0]):
        for j in range(1, matrix.shape[1]):
            h_values[0] = matrix[i - 1][j] - 1
            h_values[1] = matrix[i][j - 1] - 1
            if seq1[i - 1] == seq2[j - 1]:
                h_values[2] = matrix[i - 1][j - 1] + 1
            else:
                h_values[2] = matrix[i - 1][j - 1] - 1
            if i == d_index_i and j == d_index_j:
                matrix[i][i] = h_values[2]
            else:
                matrix[i][j] = np.max(h_values)

#@njit
def calculate_alignment_score(seq1, seq2):
    m_sum = np.zeros((len(seq1), len(seq2)))
    for i in range(0, m_sum.shape[0]):
        for j in range(0, m_sum.shape[1]):
            m = np.zeros((len(seq_t) + 1, len(seq_s) + 1))
            m[0] = -np.arange(0, m.shape[1])
            m[:, 0] = -np.arange(0, m.shape[0])
            danied_equal_fill_matrix(m, seq1, seq2, i, j)
            m_sum[i][j] = m[-1][-1]
    return m_sum


try:
    DNA_s, DNA_t = read_data()
except ValueError:
    print("\nBłąd odczytu danych!\n"
          "Możliwość znalezienia dopasowania tylko między dwoma sekwencjami! "
          "Podaj prawidłowy format danych.")
    sys.exit(0)

seq_t = str(DNA_t.seq)
seq_s = str(DNA_s.seq)

m = np.zeros((len(seq_t) + 1, len(seq_s) + 1))
m[0] = -np.arange(0, m.shape[1])
m[:, 0] = -np.arange(0, m.shape[0])
fill_matrix(m, seq_t, seq_s)

m_sum = calculate_alignment_score( seq_t, seq_s)

print(m_sum)
print(np.sum(m_sum))
print(m[-1][-1])