import sys
import numpy as np
from Bio import SeqIO
from numba import njit
from timer import Timer


def read_data(path="/home/krzysztof/Pobrane/rosalind_osym.txt"):
    """
    Odczytuje dane z pliku.
    :param path: Ścieżka do pliku.
    :return: Zwraca tablice rekordów w formacie bibilioteki biopython.
    """
    records = []
    for record in SeqIO.parse(path, "fasta"):
        records.append(record)
    return records


@njit
def fill_matrix(matrix, seq1, seq2):
    """
    Implementacja algorytmu Needlemana-Wunscha.
    :param matrix: Tablica 2D o bokach o długości - długość porównychaych sekwencji + 1.
    :param seq1: Porównywana sekwencja.
    :param seq2: Porównywana sekwencja.
    :return:
    """
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


@njit
def reverse_fill_matrix(matrix, seq1, seq2):
    """
    Implementacja algorytmu Needlemana-Wunscha. Wypełnia tablice zaczynając od ostatniej komórki (index -1 -1).
    :param matrix: Tablica 2D o bokach o długości - długość porównychaych sekwencji + 1.
    :param seq1: Porównywana sekwencja.
    :param seq2: Porównywana sekwencja.
    :return:
    """
    h_values = np.zeros((3,))
    for i in range(matrix.shape[0]-2, -1, -1):
        for j in range(matrix.shape[1]-2, -1, -1):
            h_values[0] = matrix[i + 1][j] - 1
            h_values[1] = matrix[i][j + 1] - 1
            if seq1[i] == seq2[j]:
                h_values[2] = matrix[i + 1][j + 1] + 1
            else:
                h_values[2] = matrix[i + 1][j + 1] - 1
            matrix[i][j] = np.max(h_values)


@njit
def calculate_alignment_score(alignment, reverse_alignment, seq1, seq2):
    """
    Dla każdej pary nukleotydów z sekwencji seq1 seq2 oblicza maksymalny scora, za dopasowanie, zawierające te
    dwa nukleotydy dopasowane do siebie.
    :param alignment: Macierz wypełniona zgodnie z algorytmem Needlemana Wunscha.
    :param reverse_alignment: Macierz wypełniona zgodnie z algorytmem Needlemana Wunscha, zaczynając od ostatniej
    komórki.
    :param seq1: Sekwencja nukleotydów.
    :param seq2: Sekwencja nukleotydów.
    :return: Macierz zawierającą w każdje komurce maksymalny scor, obliczony zgodnie z opisem wyżej.
    """
    m_matrix = np.zeros((len(seq1), len(seq2)))
    for i in range(0, m_matrix.shape[0]):
        for j in range(0, m_matrix.shape[1]):
            m_matrix[i][j] += 1 if seq1[i] == seq2[j] else -1
            m_matrix[i][j] += alignment[i][j]
            m_matrix[i][j] += reverse_alignment[i+1][j+1]
    return m_matrix


try:
    DNA_s, DNA_t = read_data()
except ValueError:
    print("\nBłąd odczytu danych!\n"
          "Możliwość znalezienia dopasowania tylko między dwoma sekwencjami! "
          "Podaj prawidłowy format danych.")
    sys.exit(0)

seq_t = str(DNA_t.seq)
seq_s = str(DNA_s.seq)

with Timer() as t:
    alignment = np.zeros((len(seq_t) + 1, len(seq_s) + 1))
    alignment[0] = -np.arange(0, alignment.shape[1])
    alignment[:, 0] = -np.arange(0, alignment.shape[0])
    fill_matrix(alignment, seq_t, seq_s)

    r_m = np.zeros((len(seq_t) + 1, len(seq_s) + 1))
    r_m[-1] = -np.arange(r_m.shape[1]-1, -1, -1)
    r_m[:, -1] = -np.arange(r_m.shape[0]-1, -1, -1)
    reverse_fill_matrix(r_m, seq_t, seq_s)

    m_matrix = calculate_alignment_score(alignment, r_m, seq_t, seq_s)

print('Czas rozwiązywania zadania wyniósł: %.03f sec.' % t.interval)
print(int(alignment[-1][-1]))
print(int(np.sum(m_matrix)))
