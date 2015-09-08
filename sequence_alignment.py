#!/usr/bin/env python

from matrices import Blosum62
import numpy as np

def create_affine_matrices(seq1,seq2,gap_open=-11,gap_extend=-1):
    """Creates distance matrices for sequence alignment"""
    X = np.zeros(shape=(len(seq2) + 1,len(seq1) + 1))
    Y = np.zeros(shape=(len(seq2) + 1,len(seq1) + 1))
    M = np.zeros(shape=(len(seq2) + 1,len(seq1) + 1))

    for i in range(1, len(seq2) + 1):
        X[i,0] = -float("inf")
        Y[i,0] = gap_open + (i * gap_extend)
        M[i,0] = -float("inf")
    for j in range(1, len(seq1) + 1):
        X[0,j] = gap_open + (j * gap_extend)
        Y[0,j] = -float("inf")
        M[0,j] = -float("inf")
    for j in range(1, len(seq1) + 1):
        for i in range(1, len(seq2) + 1):
            X[i,j] = max((gap_open + gap_extend + M[i,j-1]), (gap_extend + X[i,j-1]),
                    (gap_open + gap_extend + Y[i,j-1]))
            Y[i,j] = max((gap_open + gap_extend + M[i-1,j]), (gap_open + gap_extend + X[i-1,j]),
                    (gap_extend + Y[i-1,j]))
            M[i,j] = max(Blosum62(seq2[i-1],seq1[j-1]).sub_score() + M[i-1,j-1], X[i,j], Y[i,j])

    return X,Y,M

def affine_align(seq1,seq2,gap_open=-11,gap_extend=-1):
    """Performs global alignment with affine gap penalties"""
    aligned_seq1 = ''
    aligned_seq2 = ''
    X,Y,M = create_affine_matrices(seq1,seq2,gap_open,gap_extend)
    i = len(seq2)
    j = len(seq1)
    while (i > 0 or j > 0):
        if (i > 0 and j > 0 and M[i,j] == M[i-1,j-1] + Blosum62(seq2[i-1],seq1[j-1]).sub_score()):
            aligned_seq1 = seq1[j-1] + aligned_seq1
            aligned_seq2 = seq2[i-1] + aligned_seq2
            i -= 1
            j -= 1
        elif (i > 0 and M[i,j] == Y[i,j]):
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[i-1] + aligned_seq2
            i -= 1
        elif (j > 0 and M[i,j] == X[i,j]):
            aligned_seq1 = seq1[j-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            j -= 1
        else:
            break

    return aligned_seq1,aligned_seq2


if __name__ == '__main__':
    import sys
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    X,Y,M = create_affine_matrices(seq1,seq2)
    s1,s2 = affine_align(seq1,seq2)
    print s1
    print s2
