#!/usr/bin/env python

"""This module contains code related to the analysis of sequences. This
includes code related to comparisons and calculations, but also several
functions related to calculating indices and translating stand alone
nucleotide sequences into their amino acid equivalents"""

import re
import random
import sys
# Increase recursion limit for get_non_overlapping_indices function
sys.setrecursionlimit(10000)

from rediting.classes import matrices
import strings, rmath

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


def update_pos_seq_dict(seq_dict,gen,rna,codon_pos):
    """Update a dictionary based on residues and codon position"""
    if codon_pos == 1:
        search_str = 'first'
    elif codon_pos == 2:
        search_str = 'second'
    elif codon_pos == 3:
        search_str = 'third'
    for k in seq_dict.keys():
        if re.search(search_str,k):
            if re.search('gen',k):
                seq_dict[k] += gen
            else:
                seq_dict[k] += rna
    return seq_dict


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


def check_gap_percent(aa_list):
    """Checks a position for less than 50% gaps"""
    num_gaps = 0
    for aa in aa_list:
        if aa == '-':
            num_gaps += 1
    if float(num_gaps)/len(aa_list) <= 0.5:
        return True
    else:
        return False


def generate_start_sequence(length,bases):
    """Generates a new starting sequence"""
    seq = []
    for i in range(length):
        seq.append(random.choice(bases))
    return seq


def change_base(old_base,bases):
    """Ensures that new base is different from old base"""
    new_base = old_base
    while not new_base != old_base:
        new_base = random.choice(bases)
    return new_base


def mutate_sequence(start_seq,num_muts,bases):
    """This function takes a starting sequence and applies a number
    of mutations randomly along its length"""
    end = (len(start_seq) - 1)
    new_seq = start_seq[:]
    i = 0
    seen = set()
    while i < int(num_muts):
        position = random.randint(0,end)
        if position not in seen:
            new_base = change_base(start_seq[position],bases)
            new_seq[position] = new_base
            seen.add(position)
            i += 1
    return new_seq


def get_positional_information(gen_seq,rna_seq,codon_pos):
    """This function provides all relevant indices for each
    codon position in a sequence"""
    #print codon_pos
    num_edits = 0
    pos_dict = {1:[], 2:[], 3:[]}
    pos_base_dict = {
        1:[['A','G',0],['A','T',0],['A','C',0],
        ['G','A',0],['G','T',0],['G','C',0],
        ['T','A',0],['T','G',0],['T','C',0],
        ['C','A',0],['C','G',0],['C','T',0]],
        2:[['A','G',0],['A','T',0],['A','C',0],
        ['G','A',0],['G','T',0],['G','C',0],
        ['T','A',0],['T','G',0],['T','C',0],
        ['C','A',0],['C','G',0],['C','T',0]],
        3:[['A','G',0],['A','T',0],['A','C',0],
        ['G','A',0],['G','T',0],['G','C',0],
        ['T','A',0],['T','G',0],['T','C',0],
        ['C','A',0],['C','G',0],['C','T',0]]
        }
    for i,(rg,rm) in enumerate(zip(gen_seq,rna_seq)):
        if rg == '-':
            pass
        else:
            pos_dict[codon_pos].append(i)
            n = 0
            for b1,b2,x in pos_base_dict.get(codon_pos):
                if b1 == rg and b2 == rm:
                    num_edits += 1
                    pos_base_dict[codon_pos][n][2] += 1
                n += 1
            codon_pos = incr_codon_position(codon_pos)
    cumul_weights = []
    for codon_pos in pos_base_dict.keys():
        for sb,eb,w in pos_base_dict.get(codon_pos):
            cumul_weights.append(((codon_pos,sb,eb),w))
    return (num_edits,pos_dict,cumul_weights)


# Potentially obsolete
def check_weighted_base(base,base_weights):
    """Checks the chosen base by weighted random choice
    against a distribution. If 'base' is chosen returns,
    TRUE, otherwise returns FALSE"""
    base_choice = rmath.weighted_choice(base_weights)
    if base_choice == base:
        return True
    else:
        return False


def weighted_mutation(start_seq,num_muts,codon_positions,
        weights,sim_obj):
    """Supplies a list of overall weights"""
    new_seq = start_seq[:]
    i = 0
    seen = set()
    while i < int(num_muts):
        codon_pos,old_base,new_base = rmath.weighted_choice(weights)
        codon_pos = int(codon_pos)
        pos = choose_index(start_seq,old_base,codon_positions.get(codon_pos))
        if pos not in seen:
            new_seq[pos] = new_base
            sim_obj.update_transdict(old_base,new_base)
            sim_obj.update_codon_dict(codon_pos)
            seen.add(pos)
            i += 1
        else:
            pass
    return new_seq


def choose_index(seq,base,indices):
    chosen = False
    while chosen == False:
        new_index = random.sample(indices,1)[0]
        seq_base = seq[new_index]
        if seq_base == base:
            chosen = True
    return new_index


def calc_cluster_score(seq1,seq2,window_size=60):
    """Short-hand version of full edit cluster program to use
    in simulation analyses"""
    compstr1 = ''
    for i, (s1, s2) in enumerate(zip(seq1, seq2)):
        # Insertion in either
        if s1 == '-' and s2 != '-':
            compstr1 += str(0)
        elif s2 == '-' and s1 != '-':
            compstr1 += str(0)
        # No edits
        elif s1 == s2:
            compstr1 += str(0)
        # Edit
        elif s1 != s2:
            compstr1 += str(1)
        else:
            pass
    #print str(sum([int(val) for val in compstr1]))
    edit_list = []
    # Gets all full-length windows possible for length of aligned sequences
    for start,end in get_indices(compstr1, window_size):
        try:
            edit_list.append(calc_percent(compstr1, start, end, window_size))
        # This error should not get thrown, but just in case
        except(ValueError,IndexError):
            print "Error detected while adding to edit_list"
            pass
    return edit_list


def compare_aa_seqs(gen_index,gen_seq,rna_seq,ref_seq,edit_list):
    """Compares amino acids and determines whether changes result
    in increased or decreased similarity to a reference"""
    for i, (rg,rm) in enumerate(zip(gen_seq, rna_seq)):
        ident_before = 'N'
        ident_after = 'N'
        sim_before = 'N'
        sim_after = 'N'
        # Insert in mrna sequence
        if rg == '-' and rm != '-':
            pass
        # Insert in genomic sequence
        elif rg != '-' and rm == '-':
            # Don't do anything, but move the index along
            gen_index += 1
        # Gaps in both sequences, no edits possible
        elif rg == '-' and rm == '-':
            pass
        # No edits, but move the index along
        elif rg == rm:
            gen_index += 1
        # Edit detected
        elif rg != rm:
            ref_aa = ref_seq[i]
            # We only care about edits we can actually compare
            if ref_aa != '-':
                # Check for identity
                if rg == ref_aa:
                    ident_before = 'Y'
                if rm == ref_aa:
                    ident_after = 'Y'
                # Get the relevant score from Blosum62
                # Similar amino acids have positive scores
                scr_before = matrices.Blosum62(rg,ref_aa).sub_score()
                scr_after = matrices.Blosum62(rm,ref_aa).sub_score()
                # Check for similarity
                # Note that while both cannot be identical, they can
                # both be similar
                if scr_before > 0:
                    sim_before = 'Y'
                if scr_after > 0:
                    sim_after = 'Y'
                scr_diff = (scr_after - scr_before)

                edit_list.append([(gen_index+1),ref_aa,rg,rm,ident_before,
                    ident_after,sim_before,sim_after,scr_diff])
            gen_index += 1


def check_nonsynonymous_edit(cpos, gcod, mnuc):
    """Checks whether an edit on its own results in an AA change"""
    gcod_list = list(gcod)
    gcod_list[cpos-1] = mnuc
    new_gcod = "".join(gcod_list)
    old_gaa,new_gaa = translate(gcod,1), translate(new_gcod,1)
    if old_gaa != new_gaa:
        return True
    return False
