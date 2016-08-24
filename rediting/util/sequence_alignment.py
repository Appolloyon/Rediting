#!/usr/bin/env python

"""This module contains code to align two amino acid sequences along their
entire length using affine (i.e. gap open/gap extend) penalties. Part of
this code is based off code publically available at: https://github.com
/dnase/affine-gap-sequence-alignment/blob/master/alignment.py. The author
greatly appreciates code contribution from the above source."""


from rediting.classes.matrices import Blosum62
import numpy as np

def create_affine_matrices(seq1,seq2,gap_open=-11,gap_extend=-1):
    """Creates distance matrices for sequence alignment"""
    # The simplest thing to do is create zero-value matrices
    # first and then populate them afterwards, since we know
    # the shape of each matrix beforehand
    X = np.zeros(shape=(len(seq2) + 1,len(seq1) + 1))
    Y = np.zeros(shape=(len(seq2) + 1,len(seq1) + 1))
    M = np.zeros(shape=(len(seq2) + 1,len(seq1) + 1))

    for i in range(1, len(seq2) + 1):
        # "inf" is the pythonic representation of infinity
        X[i,0] = -float("inf")
        # If one sequence were all insertions
        Y[i,0] = gap_open + (i * gap_extend)
        M[i,0] = -float("inf")
    for j in range(1, len(seq1) + 1):
        # If the other sequence were all insertions
        X[0,j] = gap_open + (j * gap_extend)
        Y[0,j] = -float("inf")
        M[0,j] = -float("inf")
    # This is the main loop to populate the matrices!
    for j in range(1, len(seq1) + 1):
        for i in range(1, len(seq2) + 1):
            # Populate the first matrix
            X[i,j] = max((gap_open + gap_extend + M[i,j-1]), # Match plus gap
                    (gap_extend + X[i,j-1]), # Continue a gap
                    (gap_open + gap_extend + Y[i,j-1])) # Start a gap
            # Populate the second matrix, similar to above
            Y[i,j] = max((gap_open + gap_extend + M[i-1,j]),
                    (gap_open + gap_extend + X[i-1,j]),
                    (gap_extend + Y[i-1,j]))
            # Populate the last matrix, somewhat different
            M[i,j] = max(Blosum62(seq2[i-1],seq1[j-1]).sub_score() +\
                    M[i-1,j-1], X[i,j], Y[i,j]) # A match

    return X,Y,M

def affine_align(seq1,seq2,gap_open=-11,gap_extend=-1):
    """Performs global alignment with affine gap penalties"""
    aligned_seq1 = ''
    aligned_seq2 = ''
    # Create the matrices before traversing them
    X,Y,M = create_affine_matrices(seq1,seq2,gap_open,gap_extend)
    i = len(seq2)
    j = len(seq1)
    # "Dynamic programming" to move through each matrix
    while (i > 0 or j > 0): # Importantly, both need to hit zero!
        # If the best score is a match
        if (i > 0 and j > 0 and M[i,j] == M[i-1,j-1] +\
            Blosum62(seq2[i-1],seq1[j-1]).sub_score()):
            # Add the residue from each to the right aligned sequence
            aligned_seq1 = seq1[j-1] + aligned_seq1
            aligned_seq2 = seq2[i-1] + aligned_seq2
            # Decrement both counters
            i -= 1
            j -= 1
        # Gap in seq1
        elif (i > 0 and M[i,j] == Y[i,j]):
            # Add a gap in seq1, and a residue in seq2
            aligned_seq1 = '-' + aligned_seq1
            aligned_seq2 = seq2[i-1] + aligned_seq2
            # Only decrement counter if a residue is used up
            i -= 1
        # Gap in seq2
        elif (j > 0 and M[i,j] == X[i,j]):
            # Add a gap in seq2, and a residue in seq1
            aligned_seq1 = seq1[j-1] + aligned_seq1
            aligned_seq2 = '-' + aligned_seq2
            # Only decrement counter if a residue is use up
            j -= 1
        else:
            # If we finally hit zero on both we are done
            break

    return aligned_seq1,aligned_seq2


if __name__ == '__main__':
    # Simple test code
    import sys
    seq1 = sys.argv[1]
    seq2 = sys.argv[2]
    X,Y,M = create_affine_matrices(seq1,seq2)
    s1,s2 = affine_align(seq1,seq2)
