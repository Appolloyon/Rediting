#!/usr/bin/env python

from matrices import Blosum62
import numpy as np

def create_scoring_matrix(seq1,seq2,gap_penalty=-11):
    """Creates a 2D scoring matrix Fnm"""
    Fnm = np.zeros(shape=(len(seq1) + 1,len(seq2) + 1))
    #print Fnm
    for i in range(len(seq1) + 1):
        Fnm[i,0] = i * gap_penalty
    #print Fnm
    for j in range(len(seq2) + 1):
        Fnm[0,j] = j * gap_penalty
    #print Fnm
    for i in range(1,len(seq1) + 1):
        for j in range(1,len(seq2) + 1):
            match = Fnm[i-1,j-1] + Blosum62(seq1[i-1],seq2[j-1]).sub_score()
            delete = Fnm[i-1,j] + gap_penalty
            insert = Fnm[i,j-1] + gap_penalty
            Fnm[i,j] = max(match,delete,insert)
    #print Fnm
    return Fnm

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
        print "i is " + str(i)
        print "j is " + str(j)
        print "aligned seq1 is " + aligned_seq1
        print "aligned seq2 is " + aligned_seq2
    return (aligned_seq1,aligned_seq2)

if __name__ == '__main__':
    import sys
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    matrix_file = sys.argv[3]
    #with open(matrix_file,'w') as o:
        #o.write(create_scoring_matrix(seq1,seq2))
    np.savetxt(matrix_file,create_scoring_matrix(seq1,seq2),fmt='%d')
    s1,s2 = global_align(seq1,seq2)
    print s1
    print s2

