import sys
import os
import numpy as np
from Bio import SeqIO
from numba import njit
from timer import Timer


def read_data(path="/home/krzysztof/Pobrane/rosalind_smgb.txt"):
    records = []
    for record in SeqIO.parse(path, "fasta"):
        records.append(record)
    return records


def save_to_file(text, file_name='result'):
    path = os.path.dirname(os.path.abspath(__file__))
    itr = 0
    while os.path.exists(path + "/" + file_name + str(itr)):
        itr += 1
    path = path + "/" + file_name + str(itr)
    with open(path, "w+") as text_file:
        print(text, file=text_file)


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
    matrix[0] = np.zeros(matrix.shape[1])
    matrix[:, 0] = gap_penalty * np.arange(0, matrix.shape[0])


def set_start_h_values(h_values):
    h_values[2][0] = np.ones((h_values.shape[2]))


def find_best_index(matrix):
    index_i = np.argmax(matrix[:, -1])
    index_j = np.argmax(matrix[-1])
    if index_i < index_j:
        index_i = matrix.shape[0] - 1
    else:
        index_j = matrix.shape[1] - 1
    return index_i, index_j


@njit
def fill_begin_gap(matrix_shape, seq1, seq2, index_i, index_j):
    alignment1 = ''
    alignment2 = ''
    if index_i < matrix_shape[0] - 1:
        itr = 0
        while itr < matrix_shape[0] - index_i:
            alignment1 = seq1[itr] + alignment1
            alignment2 = '-' + alignment2
            itr += 1
    elif index_j < matrix_shape[1] - 1:
        itr = 0
        while itr < matrix_shape[1] - index_j:
            alignment2 = seq2[itr] + alignment2
            alignment1 = '-' + alignment1
            itr += 1
    return alignment1, alignment2


@njit
def find_alignment(path, index_i, index_j, seq1, seq2):
    alignment1 = ""
    alignment2 = ""
    while not index_i == index_j == 0:
        choice = path[index_i][index_j]
        if choice == 1:
            alignment1 = seq1[index_i-1] + alignment1
            alignment2 = seq2[index_j-1] + alignment2
            index_i -= 1
            index_j -= 1
        elif choice == 0:
            alignment1 = seq1[index_i-1] + alignment1
            alignment2 = '-' + alignment2
            index_i -= 1
        else:
            alignment1 = '-' + alignment1
            alignment2 = seq2[index_j-1] + alignment2
            index_j -= 1
    return alignment1, alignment2


def get_one_of_the_best(matrix, h_values, seq1, seq2):
    index_i, index_j = find_best_index(matrix)
    max_score = matrix[index_i][index_j]
    path = np.argmax(h_values, axis=0)
    alignment1, alignment2 = fill_begin_gap(matrix.shape, seq1, seq2, index_i, index_j)
    temp = find_alignment(path, index_i, index_j, seq1, seq2)
    alignment1 += temp[0]
    alignment2 += temp[1]
    return max_score, alignment1, alignment2


try:
    DNA_s, DNA_t = read_data()
except ValueError:
    print("\nBłąd odczytu danych!\n"
          "Możliwość znalezienia dopasowania tylko między dwoma sekwencjami! "
          "Podaj prawidłowy format danych.")
    sys.exit(0)

seq_t = str(DNA_t.seq)
seq_s = str(DNA_s.seq)
seq_s, seq_t = min(seq_t, seq_s), max(seq_t, seq_s)
gap_penalty = -1
matrix = np.zeros((len(seq_s) + 1, len(seq_t) + 1), dtype=np.int32)
h_values = np.zeros((3, len(seq_s) + 1, len(seq_t) + 1), dtype=np.int32)

with Timer() as t:
    set_start_matrix(matrix, gap_penalty)
    set_start_h_values(h_values)
    fill_matrix(matrix, h_values, seq_s, seq_t)
    alignment = get_one_of_the_best(matrix, h_values, seq_s, seq_t)

assert len(alignment[1]) == len(alignment[2])
print('Czas rozwiązywania zadania wyniósł: %.03f sec.' % t.interval)

text_to_file = str(alignment[0]) + '\n' + alignment[1] + "\n" + alignment[2]
save_to_file(text_to_file)
