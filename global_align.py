#!/usr/bin/env python

from matrices import Blosum62
import numpy as np

def create_scoring_matrix(seq1,seq2,gap_penalty=-11):
    """Creates a 2D scoring matrix Fnm"""
    Fnm = np.zeros(shape=(len(seq1) + 1,len(seq2) + 1))
    for i in range(len(seq1) + 1):
        Fnm[i,0] = i * gap_penalty
    print Fnm
    for j in range(len(seq2) + 1):
        Fnm[0,j] = j * gap_penalty
    print Fnm
    for i in range(1,len(seq1) + 1):
        for j in range(1,len(seq2) + 1):
            match = Fnm[i-1,j-1] + Blosum62(seq1[i-1],seq2[j-1]).sub_score()
            delete = Fnm[i-1,j] + gap_penalty
            insert = Fnm[i,j-1] + gap_penalty
            Fnm[i,j] = max(match,delete,insert)
    return Fnm

def create_affine_matrices(seq1,seq2,gap_open=-11,gap_extend=-1):
    """Creates three matrices for affine alignment"""
    M = np.zeros(shape=(len(seq1) + 1,len(seq2) + 1))
    X = np.zeros(shape=(len(seq1) + 1,len(seq2) + 1))
    Y = np.zeros(shape=(len(seq1) + 1,len(seq2) + 1))
    for i in range(1,len(seq1) + 1):
        M[i,0] = -float("inf")
        X[i,0] = -float("inf")
        Y[i,0] = gap_open + (i * gap_extend)
    for j in range(1,len(seq2) + 1):
        M[0,j] = -float("inf")
        X[0,j] = gap_open + (j * gap_extend)
        Y[0,j] = -float("inf")
    for i in range(1,len(seq1) + 1):
        for j in range(1,len(seq2) + 1):
            M[i,j] = max(Blosum62(seq1[i-1],seq2[j-1]).sub_score() +
                M[i-1,j-1], X[i,j], Y[i,j])
            X[i,j] = max(gap_open + gap_extend + M[i,j-1],
                    gap_open + X[i,j-1],
                    gap_open + gap_extend + Y[i,j-1])
            Y[i,j] = max(gap_open + gap_extend + M[i-1,j],
                    gap_open + gap_extend + X[i-1,j],
                    gap_open + Y[i-1,j])
    return (M,X,Y)
    #print M
    #print X
    #print Y


def global_align(seq1,seq2,gap_penalty=-11):
    """Performs global alignment on two amino acid sequences"""
    aligned_seq1 = ''
    aligned_seq2 = ''
    Fnm = create_scoring_matrix(seq1,seq2,gap_penalty)
    i = len(seq1) #+ 1
    j = len(seq2) #+ 1
    while (i > 0 or j > 0):
        if i > 0 and j > 0 and Fnm[i,j] == Fnm[i-1,j-1] + Blosum62(seq1[i-1],seq2[j-1]).sub_score():
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and Fnm[i,j] == Fnm[i-1,j] + gap_penalty:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        else:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
    return (aligned_seq1,aligned_seq2)

def global_affine_align(seq1,seq2,gap_open=-11,gap_extend=-1):
    """Performs global alignment with affine gap penalties on
    two amino acid sequences"""
    aligned_seq1 = ''
    aligned_seq2 = ''
    M,X,Y = create_affine_matrices(seq1,seq2,gap_open,gap_extend)
    print M
    print X
    print Y
    i = len(seq1)
    j = len(seq2)
    while (i > 0 or j > 0):
        print i
        print j
        if i > 0 and j > 0 and M[i,j] == M[i-1,j-1] + Blosum62(seq1[i-1],seq2[j-1]).sub_score():
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            i -= 1
            j -= 1
        elif i > 0 and M[i,j] == Y[i,j]:
            aligned_seq1 = seq1[i-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            i -= 1
        elif j > 0 and M[i,j] == X[i,j]:
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[j-1] + aligned_seq2
            j -= 1
        print "aligned_seq1 = " + aligned_seq1
        print "aligned_seq2 = " + aligned_seq2
    return (aligned_seq1,aligned_seq2)



if __name__ == '__main__':
    import sys
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    matrix_file = sys.argv[3]
    create_affine_matrices(seq1,seq2)
    #with open(matrix_file,'w') as o:
        #o.write(create_scoring_matrix(seq1,seq2))
    #np.savetxt(matrix_file,create_scoring_matrix(seq1,seq2),fmt='%d')
    #np.savetxt(matrix_file,create_affine_matrices(seq1,seq2))#,fmt='%d')
    #s1,s2 = global_align(seq1,seq2)
    s1,s2 = global_affine_align(seq1,seq2)
    print s1
    print s2

