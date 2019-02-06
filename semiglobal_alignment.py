import sys
import os
import numpy as np
from Bio import SeqIO
from numba import njit
from timer import Timer


def read_data(path="/home/krzysztof/Pobrane/rosalind_smgb.txt"):
    """
    Odczytuje dane z pliku.
    :param path: Ścieżka do pliku.
    :return: Zwraca tablice rekordów w formacie bibilioteki biopython.
    """
    records = []
    for record in SeqIO.parse(path, "fasta"):
        records.append(record)
    return records


def save_to_file(text, file_name='result'):
    """
    Zapisuje text do pliku file_name. Jeśli taki plik już istnieje, tworzy kolejny.
    :param text: Zmienna string, która będzie zapisana w pliku
    :param file_name: Nazwa pliku.
    :return:
    """
    path = os.path.dirname(os.path.abspath(__file__))
    itr = 0
    while os.path.exists(path + "/" + file_name + str(itr)):
        itr += 1
    path = path + "/" + file_name + str(itr)
    with open(path, "w+") as text_file:
        print(text, file=text_file)


@njit
def fill_matrix(matrix, h_values, seq1, seq2):
    """
    Implementacja algorytmu Needlemana-Wunscha.
    :param matrix: Tablica 2D o bokach o długości - długość porównychaych sekwencji + 1.
    :param h_values: Tablica 3D przechowuje wyszystkie wartości funcki H. Po trzy dla każdego wiersza tablicy.
    :param seq1: Porównywana sekwencja.
    :param seq2: Porównywana sekwencja.
    :return:
    """
    for i in range(1, matrix.shape[0]):
        for j in range(1, matrix.shape[1]):
            h_values[0][i][j] = matrix[i - 1][j] - 1
            h_values[1][i][j] = matrix[i][j - 1] - 1
            if seq1[i - 1] == seq2[j - 1]:
                h_values[2][i][j] = matrix[i - 1][j - 1] + 1
            else:
                h_values[2][i][j] = matrix[i - 1][j - 1] - 1
            matrix[i][j] = np.max(h_values[:, i, j])


def find_best_index(matrix):
    """
    Znajduje największą wartość w ostatniej kolumnie i ostatnim wierszy macierzy.
    :param matrix: Przeszukiwana macierz.
    :return: Indeksy komórki z największą wartością.
    """
    index_i = np.argmax(matrix[:, -1])
    index_j = np.argmax(matrix[-1])
    if index_i < index_j:
        index_i = matrix.shape[0] - 1
    else:
        index_j = matrix.shape[1] - 1
    return index_i, index_j


@njit
def find_alignment(path, best_i, best_j, seq1, seq2):
    """
    Odnajduje najlepsze dopasowanie, zgodnie z algorytmem Needlemana_Wunscha.
    :param path: Macierz, której komórki wskazują najlepsze kierunki poruszania się.
    :param best_i: Index, od którego zaczyna się dopasowanie.
    :param best_j: Index, od którego zaczyna się dopasowanie.
    :param seq1: Porównywana sekwencja.
    :param seq2: Porównywana sekwencja.
    :return: Dwa stringi z optymalnym dopasowaniem.
    """
    alignment1 = ""
    alignment2 = ""

    index_i, index_j = path.shape
    index_i -= 1
    index_j -= 1

    while index_i > best_i:
        alignment1 = seq1[index_i - 1] + alignment1
        alignment2 = '-' + alignment2
        index_i -= 1

    while index_j > best_j:
        alignment1 = '-' + alignment1
        alignment2 = seq2[index_j - 1] + alignment2
        index_j -= 1

    while True:
        choice = path[index_i][index_j]
        if choice == 2:
            alignment1 = seq1[index_i-1] + alignment1
            alignment2 = seq2[index_j-1] + alignment2
            index_i -= 1
            index_j -= 1
        elif choice == 1:
            alignment1 = '-' + alignment1
            alignment2 = seq2[index_j - 1] + alignment2
            index_j -= 1
        elif choice == 0:
            alignment1 = seq1[index_i - 1] + alignment1
            alignment2 = '-' + alignment2
            index_i -= 1
        if index_i == 0 and index_j == 0:
            break
    return alignment1, alignment2


def get_one_of_the_best(matrix, h_values, seq1, seq2):
    """
    Funkcja odpala funkcję find_alignment. Zawiera inicjalizacjie path, która ze względu na użycie numby nie mogła znaleść się
    w funkcji find_alignment
    :param matrix: Wypełniona macierz Needlemana-Wunscha.
    :param h_values: Wypełniona macierz z wartościami funkcji H.
    :param seq1: Porównywana sekwencja.
    :param seq2: Porównywana sekwencja.
    :return:
    """
    index_i, index_j = find_best_index(matrix)
    max_score = matrix[index_i][index_j]

    path = np.argmax(h_values, axis=0)
    path[0] = np.ones_like(path[0])
    path[:, 0] = np.zeros_like(path[:, 0])

    alignment1, alignment2 = find_alignment(path, index_i, index_j, seq1, seq2)
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
matrix = np.zeros((len(seq_s) + 1, len(seq_t) + 1), dtype=np.int32)
h_values = np.zeros((3, len(seq_s) + 1, len(seq_t) + 1), dtype=np.int32)

with Timer() as t:
    fill_matrix(matrix, h_values, seq_s, seq_t)
    alignment = get_one_of_the_best(matrix, h_values, seq_s, seq_t)

assert len(alignment[1]) == len(alignment[2])
print('Czas rozwiązywania zadania wyniósł: %.03f sec.' % t.interval)

text_to_file = str(alignment[0]) + '\n' + alignment[1] + "\n" + alignment[2]
save_to_file(text_to_file)
