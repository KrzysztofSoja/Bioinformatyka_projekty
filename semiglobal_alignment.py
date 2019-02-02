import numpy as np
from Bio import pairwise2
from Bio import SeqIO
import sys


#Wywali substitution_matrix
class NeedlemanWunsch():
    def __init__(self, sequnce1, sequnce2, substitution_matrix, penalty=-1):
        self.sequnce1 = sequnce1  # Indexowi litery w sekencji
        self.sequnce2 = sequnce2  # opowiada index o 1 większy w tablicy
        self.penalty = penalty
        self.substitution_matrix = substitution_matrix
        self.matrix = np.zeros((len(sequnce1) + 1, len(sequnce2) + 1))
        #Co to jest?
        self.directions = [[None for i in range(len(sequnce2) + 1)] \
                           for j in range(len(sequnce1) + 1)]

    #Chyba OK
    def _set_start_variables(self):
        self.matrix[0][0] = 0
        for i in range(1, self.matrix.shape[0]):
            self.matrix[i][0] = i * self.penalty
        for i in range(1, self.matrix.shape[1]):
            self.matrix[0][i] = i * self.penalty

    def _calculate_values(self, index_i, index_j):
        temp_a = self.matrix[index_i - 1][index_j] + self.penalty
        temp_b = self.matrix[index_i][index_j - 1] + self.penalty
        try:
            temp_c = self.matrix[index_i - 1][index_j - 1] \
                     + self.substitution_matrix[
                         (self.sequnce1[index_i - 1], self.sequnce2[index_j - 1])]
        except KeyError:
            temp_c = -np.inf
        return [(temp_a, (-1, 0)), (temp_b, (0, -1)), (temp_c, (-1, -1))]

    def _find_the_best_fit(self, values):
        answer = [max(values)]
        values.remove(answer[0])
        for i in values[:-1]:
            if i[0] == answer[0][0]:
                answer.append(i)
        return answer

    def start(self):
        self._set_start_variables()
        for i in range(1, self.matrix.shape[0]):
            for j in range(1, self.matrix.shape[1]):
                temp = self._calculate_values(i, j)
                temp = self._find_the_best_fit(temp)
                self.matrix[i, j] = temp[0][0]
                self.directions[i][j] = [d[1] for d in temp]

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


print(DNA_s.seq)
print(DNA_t.seq)

NeedlemanWunsch.get_fit()

#print(pairwise2.align.globalxx(human_str, rat_str))