#!/usr/bin/env python

import math

def gulp(string, start, gulp_size):
    """get substrings of a string"""
    gulpstr = ''
    chars = string[start:start+gulp_size]
    for char in chars:
        gulpstr += char
    return gulpstr

def compare_seqs(seq1, seq2, num_equal):
    """compare substrings to determine start of alignment"""
    equal = 0
    for i, (r1, r2) in enumerate(zip(seq1, seq2)):
        if i == 0:  #terminal residue
            if r1 != '-' and r2 != '-':  #neither should be a gap
                if r1 == r2:
                    equal += 1
                else:
                    pass
            else:
                return False
        else:  #other residues
            if r1 == r2:
                equal += 1
            else:
                pass
    if equal >= num_equal:  #arbitrary threshold
        return True
    else:
        return False

def nonblank_lines(f):
    """skip blank lines"""
    for l in f:
        line = l.strip('\n')
        if line:
            yield line

def sanitize(seq):
    """remove gap characters"""
    nseq = ''
    for char in seq:
        if char == '-':
            pass
        else:
            nseq += char
    return nseq

def calc_gc(string):
    """calculates GC content of a string"""
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

def build_seqdict(infile, seqdict):
    """Builds a dictionary of header/sequence key value pairs"""
    with open(infile,'U') as f:
        for line in nonblank_lines(f):
            line = line.strip('\n')
            if line.startswith(">"):
                line = line.strip(">")
                ID = line
                seqdict[ID] = ''
            else:
                seqdict[ID] += line
    return seqdict

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
    """Returns a list of start and end coordinates for windows"""
    indices = []
    w = int(window_size)
    for i in range(len(string) - w):
        try:
            index_low = i
            index_high = i + w
            indices.append([index_low, index_high])
        except(ValueError,IndexError):
            pass
    return indices

def get_non_overlapping_indices(string, start, size, indices):
    """Returns non-overlapping start and end coordinates"""
    size = int(size)
    end = start + size
    if end > len(string):
        return indices
    else:
        indices.append([start,end])
        start = end
        get_non_overlapping_indices(string, start, size, indices)
    return indices

def calc_mean(values):
    """Calculates mean of a set of values"""
    sum = 0.0
    for value in values:
        value = float(value)
        sum += value
    mean = sum/(float(len(values)))
    return mean

def calc_pearson(xvalues, yvalues, xmean, ymean):
    """Calculates a pearson correlation value"""
    N = len(xvalues)
    num = 0.0
    xdenom = 0.0
    ydenom = 0.0
    for i in range(N):
        num += ((xvalues[i] - xmean) * (yvalues[i] - ymean))
        xdenom += ((xvalues[i] - xmean)**2)
        ydenom += ((yvalues[i] - ymean)**2)
    denom = (math.sqrt(xdenom)) * (math.sqrt(ydenom))
    return num/denom

def calc_tvalue(PC, N):
    """Calculates a t value for a given pearson coefficient"""
    return abs((PC * math.sqrt(N-2))/(math.sqrt(1-(PC**2))))

def polyT(string):
    """find stretches of 4 or more T's in a row"""
    i = 0
    while i <= len(string) - 4:
        polyt = gulp(string, i, 4)
        if polyt == 'TTTT':
            return True
        else:
            pass
        i += 1

def polyTpercent(string, percent):
    """find stretches of X% T"""
    tcounter = 0
    for char in string:
        if char == 'T':
            tcounter += 1
    if tcounter >= (percent/10):
        return True
    else:
        pass

def ispolyTpercent(plist, percent):
    """check list elements for at least one polyT stretch"""
    for e in plist:
        if polyTpercent(e, percent):
            return True
        else:
            pass
    return False

def incr_codon_position(codon_pos):
    """increments codon counter"""
    if codon_pos < 3:
        codon_pos += 1
    else:
        codon_pos = 1
    return codon_pos

def calculate_codons(nuc_seq,codon_pos):
    """returns a list of codons based on reading frame"""
    codon_list = []
    # shift the reading frame based on codon position
    if int(codon_pos) == 1:
        test_seq = nuc_seq
    elif int(codon_pos) == 2:
        test_seq = nuc_seq[2:]
    elif int(codon_pos) == 3:
        test_seq = nuc_seq[1:]
    # split remaining sequence into codons
    for start,end in get_non_overlapping_indices(test_seq,0,3,indices=[]):
        codon_list.append(test_seq[start:end])
    return codon_list

def translate(nuc_seq,codon_pos):
    """returns the amino acid translation of a nucleotide sequence"""
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
    # special case when sequence starts with an insertion
    if nuc_seq[0] == '-':
        codon_pos = incr_codon_position(codon_pos)
    nuc_seq = sanitize(nuc_seq) # remove all gap characters prior to translation
    for codon in calculate_codons(nuc_seq,codon_pos):
        codon_str += codon + ', '
        aa = '-'
        for k in aa_dict.keys():
            for e in aa_dict.get(k):
                if e == codon:
                    aa = k
        if aa == 'STOP':
            aa = '-'
        aa_seq += aa
    return aa_seq

