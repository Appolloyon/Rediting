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

def calc_mean(values):
    """Calculates mean of a set of values"""
    sum = 0.0
    for value in values:
        value = float(value)
        sum += value
    mean = sum/(float(len(values)))
    return mean

def calc_pearson(xvalues, yvalues, xmean, ymean):
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
    return abs((PC * math.sqrt(N-2))/(math.sqrt(1-(PC**2))))

