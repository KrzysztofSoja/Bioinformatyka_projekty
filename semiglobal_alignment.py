import numpy as np
from Bio import pairwise2
from Bio import SeqIO
import sys
from tqdm import tqdm
from numba import jit, njit


class NeedlemanWunsch():

    directions = {0: (-1, 0),
                  1: (-1, -1),
                  2: (0, -1)}

    def __init__(self, sequnce1, sequnce2, penalty=-1):
        self.sequnce1 = sequnce1
        self.sequnce2 = sequnce2
        self.penalty = penalty
        self.matrix = np.zeros((len(sequnce1) + 1, len(sequnce2) + 1), dtype=int)
        self.h_values = np.zeros((3, len(sequnce1) + 1, len(sequnce2) + 1), dtype=int)

    #https://numba.pydata.org/numba-doc/latest/reference/pysupported.html
    #Koniecznie przeczytaj
    @njit
    def _set_start_variables(self):
        self.matrix[0] = self.penalty * np.arange(0, self.matrix.shape[1])
        self.matrix[:, 0] = self.penalty * np.arange(0, self.matrix.shape[0])


    def _set_the_best_values(self, index_i, index_j):
        pass


    def _find_the_best_fit(self, values):
        answer = [max(values)]
        values.remove(answer[0])
        for i in values[:-1]:
            if i[0] == answer[0][0]:
                answer.append(i)
        return answer

    @jit(nonpython=True)
    def fill_matrix(self):
        self._set_start_variables()
        for i in tqdm(range(1, self.matrix.shape[0])):
            for j in range(1, self.matrix.shape[1]):
                self.h_values[0][i][j] = self.matrix[i - 1][j] + self.penalty
                self.h_values[1][i][j] = self.matrix[i][j - 1] + self.penalty
                if self.sequnce1[i - 1] == self.sequnce2[j - 1]:
                    self.h_values[2][i][j] = self.matrix[i - 1][j - 1] + 1
                else:
                    self.h_values[2][i][j] = self.matrix[i - 1][j - 1] - 1
                self.matrix[i][j] = max(self.h_values[:, i, j])

    #Co tu się dzieje???
    def _get_fit(self, index_i, index_j, number_of_fit, all_fit):
        vectors = self.directions[index_i][index_j]
        for i, vector in enumerate(vectors):
            new_index_i = index_i + vector[0]
            new_index_j = index_j + vector[1]
            if i > 0:
                all_fit.insert(number_of_fit + i, list(all_fit[number_of_fit]))
            if new_index_i >= 1 and new_index_j >= 1:
                all_fit[number_of_fit + i].append(self._get_char_to_fit((
                    new_index_i, new_index_j), vector))
                all_fit = self._get_fit(new_index_i, new_index_j,
                                        number_of_fit + i, all_fit)
        return all_fit

    #???
    def _get_char_to_fit(self, index, vector):
        if vector == (-1, -1):
            return (self.sequnce1[index[0] - 1], self.sequnce2[index[1] - 1])
        elif vector == (-1, 0):
            return (self.sequnce1[index[0] - 1], '-')
        elif vector == (0, -1):
            return ('-', self.sequnce2[index[1] - 1])

    def get_fit(self):
        all_fit = [[]]
        index_i = len(self.directions) - 1
        index_j = len(self.directions[1]) - 1
        all_fit[0].append((self.sequnce1[index_i - 1],
                           self.sequnce2[index_j - 1]))
        all_fit = self._get_fit(index_i, index_j, 0, all_fit)

        for fit in all_fit:
            fit.reverse()

        return all_fit

    def print_as_string(tab_of_tuple):
        str_seq1 = ''
        str_seq2 = ''
        for i in tab_of_tuple:
            str_seq1 = str_seq1 + i[0]
            str_seq2 = str_seq2 + i[1]
        print(str_seq1)
        print(str_seq2)


def read_data(path="/home/krzysztof/Pobrane/rosalind_smgb.txt"):
    records = []
    for record in SeqIO.parse(path, "fasta"):
        records.append(record)
    return records

try:
    DNA_s, DNA_t = read_data()
except ValueError:
    print("\nBłąd odczytu danych!\n"
          "Możliwość znalezienia dopasowania tylko między dwoma sekwencjami! "
          "Podaj prawidłowy format danych.")
    sys.exit(0)


@jit(nonpython=True)
def fill_matrix(matrix, h_values, seq1, seq2):
    for i in tqdm(range(1, matrix.shape[0])):
        for j in range(1, matrix.shape[1]):
            h_values[0][i][j] = matrix[i - 1][j] - 1
            h_values[1][i][j] = matrix[i][j - 1] - 1
            if seq1[i - 1] == seq2[j - 1]:
                h_values[2][i][j] = matrix[i - 1][j - 1] + 1
            else:
                h_values[2][i][j] = matrix[i - 1][j - 1] - 1
            matrix[i][j] = max(h_values[:, i, j])


print(DNA_s.seq)
print(DNA_t.seq)

test = NeedlemanWunsch(DNA_s.seq, DNA_t.seq)
print(test.matrix.shape)
test._set_start_variables()
fill_matrix(test.matrix, test.h_values, test.sequnce1, test.sequnce2)
print(test.matrix)

#print(pairwise2.align.globalxx(human_str, rat_str))