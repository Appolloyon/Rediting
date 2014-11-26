#!/usr/bin/env python

"""
This program reads in aligned files with mRNA and DNA sequences of a gene from
one organism and the corresponding genome sequence from Emiliania. It calculates
the z score of editing between the transcript and the DNA sequence, and then
compares the sequence identity between the two genome sequences over the same
range. Both of these values are graphed.

Changelog:
Created August 18, 2014
Last updated September 5, 2014
    -updated program comments
    -added changelog section
"""

import os
import sys
import re
import matplotlib.pyplot as plt

window_size = float(sys.argv[1])
file_list = sys.argv[2:]

def nonblank_lines(f):
    for l in f:
        line = l.strip('\n')
        if line:
            yield line

def calc_zscore(string, start, end, window_size):
    chars = string[start:end]
    sum = 0.0
    w = window_size
    for char in chars:
        sum += float(char)
    zscore = float((sum/w)*100)
    return zscore

def calc_ident(string, start, end, window_size):
    chars = string[start:end]
    sum = 0.0
    w = window_size
    for char in chars:
        sum += float(char)
    pident = float((sum/w)*100)
    return pident

def get_indices(string, window_size):
    indices = []
    w = int(window_size)
    for i in range(len(string)):
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

def calc_mean(zlist):
    sum = 0.0
    for d in zlist:
        d = float(d)
        sum += d
    mean = sum/(float(len(zlist)))
    return mean

def calc_stddev(zlist, mean):
    sum = 0.0
    m = mean
    for z in zlist:
        v = float(z)
        v2 = ((v - m)**2)
        sum += v2
    stddev = ((sum/float(len(zlist)))**0.5)
    return stddev

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

        compstr1 = ''
        for i, (res1, res2) in enumerate(zip(newseq1, newseq2)):
            if (res1 == '-' or res2 == '-') or res1 == res2:
                compstr1 += str(0)
            elif res1 != res2:
                compstr1 += str(1)
            else:
                pass

        zlist = []
        for s,e in get_indices(compstr1, window_size):
            try:
                zlist.append(calc_zscore(compstr1, s, e, window_size))
            except(ValueError,IndexError):
                pass

        compstr2 = ''
        for i, (res1, res2) in enumerate(zip(newseq2, newseq3)):
            if res1 == res2:
                compstr2 += str(1)
            elif res1 != res2:
                compstr2 += str(0)
            else:
                pass

        ilist = []
        for s,e in get_indices(compstr2, window_size):
            try:
                ilist.append(calc_ident(compstr2, s, e, window_size))
            except(ValueError,IndexError):
                pass

        mean = calc_mean(zlist)
        stddev2 = 2 * (calc_stddev(zlist, mean))
        stddev3 = 3 * (calc_stddev(zlist, mean))
        #print stddev
        stddevlist2 = []
        stddevlist3 = []

        for i in range(len(zlist)):
            stddevlist2.append(stddev2)
            stddevlist3.append(stddev3)

        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot([z for z in zlist], color='b')
        ax1.plot([stddev for stddev in stddevlist2], color='k')
        ax1.plot([stddev for stddev in stddevlist3], color='k')
        ax2.plot([i for i in ilist], color='m')
        plt.title('%s' % (name))
        ax1.set_xlabel('nucleotide')
        ax1.set_ylabel('z score', color='b')
        ax2.set_ylabel('sequence identity', color='m')
        plt.savefig('%s.pdf' % (name))
        plt.close()
