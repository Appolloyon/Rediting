#!/usr/bin/env python

"""
This program reads in aligned files with mRNA and DNA sequences of a gene from
one organism and the corresponding genome sequence from Emiliania. It calculates
the percent edits between the transcript and the DNA sequence, and then
compares the sequence identity between the two genome sequences over the same
range. Both of these values are graphed and a pearson's correlation is calculated.

Changelog:
Created June 25, 2015
"""

import os
import sys
import re
import math
import matplotlib.pyplot as plt

window_size = float(sys.argv[1])
file_list = sys.argv[2:]

def nonblank_lines(f):
    for l in f:
        line = l.strip('\n')
        if line:
            yield line

def compare_seqs(seq1, seq2):
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
    if equal >= 7:  #arbitrary threshold
        return True
    else:
        return False

def calc_percent(string, start, end, window_size):
    chars = string[start:end]
    sum = 0.0
    w = window_size
    for char in chars:
        sum += float(char)
    percent = float((sum/w)*100)
    return percent

def get_indices(string, window_size):
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

def gulp(string, start, gulp_size):
    gulpstr = ''
    chars = string[start:start+gulp_size]
    for char in chars:
        gulpstr += char
    return gulpstr

def calc_mean(values):
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


for file in file_list:
    name = ((os.path.basename((file)).strip(".afa")))
    with open(file,'U') as f:
        seqdict={}
        for curline in nonblank_lines(f):
            curline = curline.strip('\n')
            if curline.startswith(">"):
                curline = curline.strip(">")
                ID = curline
                seqdict[ID] = ''
            else:
                seqdict[ID] += curline

        for k in seqdict.keys():
            if re.search('mRNA',k):
                seq1 = seqdict.get(k)
            elif re.search('Emiliania',k):
                seq3 = seqdict.get(k)
            else:
                seq2 = seqdict.get(k)

        i = 0
        while not compare_seqs((gulp(seq1, i, 9)), (gulp(seq2, i, 9))):  #start of alignment
            i += 1
        print i
        j = 0
        while not compare_seqs((gulp(seq1[::-1], j, 9)), (gulp(seq2[::-1], j, 9))):  #end of alignment
            j += 1
        print j

        newseq1 = seq1[i:(len(seq1)-j)]
        newseq2 = seq2[i:(len(seq2)-j)]
        newseq3 = seq3[i:(len(seq3)-j)]

        # calculates % edits
        compstr1 = ''
        for i, (res1, res2) in enumerate(zip(newseq1, newseq2)):
            if (res1 == '-' or res2 == '-') or res1 == res2:
                compstr1 += str(0)
            elif res1 != res2:
                compstr1 += str(1)
            else:
                pass

        edit_list = []
        for s,e in get_indices(compstr1, window_size):
            #print s
            #print e
            #print compstr1[s:e]
            try:
                edit_list.append(calc_percent(compstr1, s, e, window_size))
            except(ValueError,IndexError):
                pass

        # calculates % sequence identity
        compstr2 = ''
        for i, (res1, res2) in enumerate(zip(newseq2, newseq3)):
            if res1 == res2:
                compstr2 += str(1)
            elif res1 != res2:
                compstr2 += str(0)
            else:
                pass

        identity_list = []
        for s,e in get_indices(compstr2, window_size):
            try:
                identity_list.append(calc_percent(compstr2, s, e, window_size))
            except(ValueError,IndexError):
                pass

        edit_mean = calc_mean(edit_list)
        average_list = []
        for i in range(len(edit_list)):
            average_list.append(edit_mean)
        identity_mean = calc_mean(identity_list)

        edits_above_average = 0.0
        total_edits = 0.0
        for edit in edit_list:
            if edit > edit_mean:
                total_edits += edit
                edits_above_average += edit
            else:
                total_edits += edit

        percent_above_average_edits = (edits_above_average/total_edits) * 100
        print percent_above_average_edits

        PC = calc_pearson(edit_list, identity_list, edit_mean, identity_mean)
        print PC
        print calc_tvalue(PC, len(edit_list))

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot([e for e in edit_list], color='b')
        ax1.plot([mean for mean in average_list], color='k')
        ax2.plot([i for i in identity_list], color='m')
        plt.title('%s' % (name))
        ax1.set_xlabel('sliding window position')
        ax1.set_ylabel('% edits', color='b')
        ax2.set_ylabel('sequence identity', color='m')
        plt.savefig('%s.pdf' % (name))
        plt.close()
