from Bio import SeqIO
import numpy as np
from numba import njit
from timer import Timer


def read_data(path="/home/krzysztof/Pobrane/rosalind_mgap.txt"):
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
def fill_matrix(score, gaps, seq1, seq2):
    """
    Wypełnia dwie macierze: score - zgodnie z algorytmem Needlemana-Wunscha oraz macierz gaps, której każda komórka
    zawiera liczbę luk, które pojawieją się od punktu 0 0 do danej komórki macierzy.
    :param score: Tablica 2D o bokach o długości - długość porównychaych sekwencji + 1. Z pierwszą kolumą oraz pierwszym wierszem
    zainicjalizowanym odpowiednią wartością kary za przerwę.
    :param gaps: Tablica 2D o bokach o długości - długość porównychaych sekwencji + 1.
    :param seq1: Porównywana sekwencja.
    :param seq2: Porównywana sekwencja.
    :return:
    """
    h_values = np.zeros((3,))
    for i in range(1, score.shape[0]):
        for j in range(1, score.shape[1]):
            h_values[0] = score[i - 1][j] - 1
            h_values[1] = score[i][j - 1] - 1
            if seq1[i - 1] == seq2[j - 1]:
                h_values[2] = score[i - 1][j - 1] + 1
            else:
                h_values[2] = -np.inf
            score[i][j] = np.max(h_values)
            arrow = np.argmax(h_values)
            if arrow == 0:
                gaps[i][j] = gaps[i - 1][j] + 1
            elif arrow == 1:
                gaps[i][j] = gaps[i][j - 1] + 1
            elif arrow == 2:
                gaps[i][j] = gaps[i - 1][j - 1]


with Timer() as t:
    seq1, seq2 = read_data()
    seq1, seq2 = str(seq1.seq), str(seq2.seq)

    alignment_matrix = np.zeros((len(seq1) + 1, len(seq2) + 1))
    alignment_matrix[0] = -np.arange(0, alignment_matrix.shape[1])
    alignment_matrix[:, 0] = -np.arange(0, alignment_matrix.shape[0])

    gaps = np.zeros((len(seq1) + 1, len(seq2) + 1))
    gaps[0] = np.arange(0, gaps.shape[1])
    gaps[:, 0] = np.arange(0, gaps.shape[0])

    fill_matrix(alignment_matrix, gaps, seq1, seq2)

print("Zadanie wykonane w czasie: %.3f sekund" % t.interval)
print(int(gaps[-1][-1]))