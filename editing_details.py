#!/usr/bin/env python

"""
This program is a companion to the 'editing_count.py' series of scripts.
Instead of graphing values, this program calculates some statistics about the
editing events noted by the graphs. Between the DNA and mRNA sequence, the
frequency of different base transitions is noted, i.e. A to T, G to C, etc,
and the orientation of the change is from DNA to mRNA.
"""

import os
import sys
import re

file_list = sys.argv[1:]

def nonblank_lines(f):
    for l in f:
        line = l.strip('\n')
        if line:
            yield line

def gulp(string, start, gulp_size):
    gulpstr = ''
    chars = string[start:start+gulp_size]
    for char in chars:
        gulpstr += char
    return gulpstr

def base_transition(seq1, seq2):
    transdict = {'a_t':0, 'a_g':0, 'a_c':0,
            't_a':0, 't_g':0, 't_c':0,
            'g_a':0, 'g_t':0, 'g_c':0,
            'c_a':0, 'c_t':0, 'c_g':0}
    for i, (b1, b2) in enumerate(zip(seq1, seq2)):
        if b1 == b2:
            pass
        elif b1 == "A" and b2 == "T":
            transdict['a_t'] += 1
        elif b1 == "A" and b2 == "G":
            transdict['a_g'] += 1
        elif b1 == "A" and b2 == "C":
            transdict['a_c'] += 1
        elif b1 == "T" and b2 == "A":
            transdict['t_a'] += 1
        elif b1 == "T" and b2 == "G":
            transdict['t_g'] += 1
        elif b1 == "T" and b2 == "C":
            transdict['t_c'] += 1
        elif b1 == "G" and b2 == "A":
            transdict['g_a'] += 1
        elif b1 == "G" and b2 == "T":
            transdict['g_t'] += 1
        elif b1 == "G" and b2 == "C":
            transdict['g_c'] += 1
        elif b1 == "C" and b2 == "A":
            transdict['c_a'] += 1
        elif b1 == "C" and b2 == "T":
            transdict['c_t'] += 1
        elif b1 == "C" and b2 == "G":
            transdict['c_g'] += 1
        else:
            pass
    return transdict

def calc_gc(string):
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

for file in file_list:
    name = ((os.path.basename((file)).strip(".afa")))
    outname = name + "_editing_stats.txt"
    with open(file,'U') as f, open(outname, 'w') as o:
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
        while gulp(seq1, i, 3) != gulp(seq2, i, 3):
            i += 1
        #print i

        j = 0
        while gulp(seq1[::-1], j, 3) != gulp(seq2[::-1], j, 3):
            j += 1
        #print i
        newseq1 = seq1[i:(len(seq1)-j)]
        newseq2 = seq2[i:(len(seq2)-j)]
        newseq3 = seq3[i:(len(seq3)-j)]

        edited_res = 0
        for i, (res1, res2) in enumerate(zip(newseq1, newseq2)):
            if (res1 == '-' or res2 == '-') or res1 == res2:
                edited_res += 0
            elif res1 != res2:
                edited_res += 1
            else:
                pass

        transdict = base_transition(newseq2, newseq1)
        DNA_gc = calc_gc(newseq2)
        mRNA_gc = calc_gc(newseq1)

        o.write("Total number of unequal residues: {}\n".format(edited_res))
        o.write("Number of transitions: \n")
        o.write("A to T: {}\n".format(transdict.get('a_t')))
        o.write("A to G: {}\n".format(transdict.get('a_g')))
        o.write("A to C: {}\n".format(transdict.get('a_c')))
        o.write("T to A: {}\n".format(transdict.get('t_a')))
        o.write("T to G: {}\n".format(transdict.get('t_g')))
        o.write("T to C: {}\n".format(transdict.get('t_c')))
        o.write("G to A: {}\n".format(transdict.get('g_a')))
        o.write("G to T: {}\n".format(transdict.get('g_t')))
        o.write("G to C: {}\n".format(transdict.get('g_c')))
        o.write("C to A: {}\n".format(transdict.get('c_a')))
        o.write("C to T: {}\n".format(transdict.get('c_t')))
        o.write("C to G: {}\n".format(transdict.get('c_g')))
        o.write("GC content before editing: {}\n".format(DNA_gc))
        o.write("GC content after editing: {}\n".format(mRNA_gc))
