#!/usr/bin/env python

"""This module contains code related to the analysis of sequences. This
includes code related to comparisons and calculations, but also several
functions related to calculating indices and translating stand alone
nucleotide sequences into their amino acid equivalents"""


import strings

def compare_seqs(seq1, seq2, num_equal):
    """Compare substrings to determine start of alignment"""
    equal = 0
    for i, (r1, r2) in enumerate(zip(seq1, seq2)):
        if i == 0:  # Terminal residue
            if r1 != '-' and r2 != '-':  # Neither should be a gap
                if r1 == r2:
                    equal += 1
                else:
                    pass
            else:
                # If the terminal residue includes a gap then we
                # return False immediately
                return False
        else:  # Non-terminal residues
            if r1 == r2:
                equal += 1
            else:
                pass
    if equal >= num_equal:  # Arbitrary threshold
        # True indicates that we have sufficient similarity to
        # count the sequences as well aligned
        return True
    else:
        return False


def calc_gc(string):
    """Calculates GC content of a string"""
    GC = 0
    AT = 0
    for char in string:
        if char == "G" or char == "C":
            GC += 1
        elif char == "A" or char == "T":
            AT += 1
        else:
            pass
    gc_content = (GC/float(GC + AT)) * 100
    return gc_content


def calc_percent(string, start, end, window_size):
    """Calculates percent value of a given string"""
    chars = string[start:end]
    sum = 0.0
    w = window_size
    for char in chars:
        sum += float(char)
    percent = float((sum/w)*100)
    return percent


def get_indices(string, window_size):
    """Returns a list of start and end coordinates for overlapping windows"""
    indices = []
    w = int(window_size)
    # Range takes up to the last number, therefore use 'w-1'
    # This ensures that we don't add windows at the end of the sequence
    # that are less than the window length
    for i in range(len(string) - (w-1)):
        try:
            index_low = i
            index_high = i + w
            indices.append([index_low, index_high])
        # We shouldn't throw this error, but just in case
        except(ValueError,IndexError):
            pass
    return indices


def get_non_overlapping_indices(string, start, size, indices):
    """Returns non-overlapping start and end coordinates"""
    size = int(size)
    end = start + size
    if end > len(string):
        # Base case for the recursion!
        # This ensures that no windows are shorter
        # than the specified window length
        return indices
    else:
        indices.append([start,end])
        # Hence, "non-overlapping"
        start = end
        # Recurse
        get_non_overlapping_indices(string, start, size, indices)
    return indices


def polyT(string):
    """Find stretches of 4 or more T's in a row"""
    i = 0
    while i <= len(string) - 4:
        # Take four bp snapshots
        polyt = strings.gulp(string, i, 4)
        if polyt == 'TTTT':
            # If even one stretch is true, evaluates True
            return True
        else:
            pass
        i += 1


def polyTpercent(string, percent):
    """Find stretches of X% T"""
    tcounter = 0
    for char in string:
        if char == 'T':
            tcounter += 1
    if tcounter >= (percent/10):
        return True
    else:
        pass


def incr_codon_position(codon_pos):
    """Increments codon counter"""
    # Implementation identical to that used in classes
    if codon_pos < 3:
        codon_pos += 1
    else:
        codon_pos = 1
    return codon_pos


def calculate_codons(nuc_seq,codon_pos):
    """Returns a list of codons based on reading frame"""
    codon_list = []
    # Shift the sequence into RF +1 based on codon position
    if int(codon_pos) == 1:
        test_seq = nuc_seq
    elif int(codon_pos) == 2:
        test_seq = nuc_seq[2:]
    elif int(codon_pos) == 3:
        test_seq = nuc_seq[1:]
    # Split remaining sequence into codons
    # We have to use non-overlapping indices here
    for start,end in get_non_overlapping_indices(test_seq,0,3,indices=[]):
        codon_list.append(test_seq[start:end])
    return codon_list


def translate(nuc_seq,codon_pos):
    """Returns the amino acid translation of a nucleotide sequence"""
    aa_dict = {
            'F':{'TTT','TTC'},
            'L':{'TTA','TTG','CTT','CTC','CTA','CTG'},
            'I':{'ATT','ATC','ATA'},
            'M':{'ATG'},
            'V':{'GTT','GTC','GTA','GTG'},
            'S':{'TCT','TCC','TCA','TCG','AGT','AGC'},
            'P':{'CCT','CCC','CCA','CCG'},
            'T':{'ACT','ACC','ACA','ACG'},
            'A':{'GCT','GCC','GCA','GCG'},
            'Y':{'TAT','TAC'},
            'H':{'CAT','CAC'},
            'Q':{'CAA','CAG'},
            'N':{'AAT','AAC'},
            'K':{'AAA','AAG'},
            'D':{'GAT','GAC'},
            'E':{'GAA','GAG'},
            'C':{'TGT','TGC'},
            'W':{'TGG'},
            'R':{'CGT','CGC','CGA','CGG','AGA','AGG'},
            'G':{'GGT','GGC','GGA','GGG'},
            'STOP':{'TAA','TAG','TGA'}
            }
    aa_seq = ''
    codon_str = ''
    # Remove all gap characters prior to translation
    nuc_seq = strings.sanitize(nuc_seq)
    for codon in calculate_codons(nuc_seq,codon_pos):
        # Break up into codons
        codon_str += codon + ', '
        aa = '-'
        # Translate all codons
        for k in aa_dict.keys():
            for e in aa_dict.get(k):
                if e == codon:
                    aa = k
        # STOP codons equivalent to gaps
        if aa == 'STOP':
            aa = '-'
        aa_seq += aa
    return aa_seq

def check_indices(start,stop):
    """Returns True if indices are sensible and False if not"""
    # This shouldn't happen, but just in case
    if stop < start:
        return False
    # If the gap is less than a codon, return False
    elif (stop - start + 1) < 3:
        return False
    # If the gap is not divisible by three return False
    elif (stop - start + 1) % 3 != 0:
        return False
    else:
        # Otherwise the indices should be fine
        return True

def expand_indices(start,stop,indices):
    """Returns a single list of index values
    for multiple start, stop pairs"""
    # Since we go up to stop, have to increment by 1
    stop = stop + 1
    while start < stop:
        # Appends all values between, and including, the
        # start and stop values to a single list
        indices.append(start)
        start += 1

def close_gap(i,codon_pos):
    """Returns an index for the end of a gap"""
    # Depending on where a gap ends, we might need to
    # move the index back to keep the right reading frame
    if codon_pos == 1:
        end = i - 1
    elif codon_pos == 2:
        end = i - 2
    elif codon_pos == 3:
        end = i - 3
    return end

def trim_sequence(seq,indices):
    """Returns a sequence lacking residues corresponding to indices"""
    new_seq = ''
    # Use enumerate to match index to each base
    for i,res in enumerate(seq):
        if i in indices:
            # Just don't add it, same as deleting
            pass
        else:
            new_seq += res
    return new_seq
