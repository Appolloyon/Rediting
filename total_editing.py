#!/usr/bin/env python

"""
Changelog
---------
Author: Christen Klinger
Created: November 25, 2014
Last Updated: November 25, 2014
"""

import re
import argparse
import os

parser = argparse.ArgumentParser(
    description = """Calculates editing stats between genomic/RNA sequences""",
    epilog = """This program assumes that the genomic and RNA sequence for a
    given gene are provided in the same file in aligned FASTA format (.afa). It
    also assumes that the files have the string 'gen' or 'mrna' in the filename,
    for genomic and mRNA sequence, respectively. Output is given as a series of
    csv values for various pertinent stats.""")
parser.add_argument('infiles', nargs='+', help='list of aligned infiles')
args = parser.parse_args()

def gulp(string, start, gulp_size):
    gulpstr = ''
    chars = string[start:start+gulp_size]
    for char in chars:
        gulpstr += char
    return gulpstr

def compare_seqs(seq1, seq2):
    equal = 0
    for i, (r1, r2) in enumerate(zip(seq1, seq2)):
        if i == 0:
            if r1 != '-' and r2 != '-':
                if r1 == r2:
                    equal += 1
                else:
                    pass
            else:
                return False
        else:
            if r1 == r2:
                equal += 1
            else:
                pass
    if equal >= 5:
        return True
    else:
        return False

def nonblank_lines(f):
    for l in f:
        line = l.strip('\n')
        if line:
            yield line


for file in args.infiles:
    name = (os.path.basename(file))
    with open(file,'U') as f:
        seqdict={}
        for line in nonblank_lines(f):
            line = line.strip('\n')
            if line.startswith(">"):
                line = line.strip(">")
                ID = line
                seqdict[ID] = ''
            else:
                seqdict[ID] += line

        for k in seqdict.keys():
            if re.search('mrna',k):
                mseq = seqdict.get(k)
            else:
                gseq = seqdict.get(k)

    i = 0
    while not compare_seqs((gulp(mseq, i, 6)), (gulp(gseq, i, 6))):
        i += 1
    j = 0
    while not compare_seqs((gulp(mseq[::-1], j, 6)), (gulp(gseq[::-1], j, 6))):
        j += 1
    newmseq = mseq[i:(len(mseq)-j)]
    newgseq = gseq[i:(len(gseq)-j)]
    print newmseq
    print newgseq
