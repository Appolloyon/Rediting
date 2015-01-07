#!/usr/bin/env python

"""
Changelog
---------
Author: Christen Klinger
Created: January 6, 2015
Last Updated: January 6, 2015
"""

import re
import argparse
from classes import CodonPair

parser = argparse.ArgumentParser(
    description = """Calculates codon preference between genomic/RNA sequences""",
    epilog = """This program assumes that the genomic and RNA sequence for a
    given gene are provided in the same file in aligned FASTA format (.afa). It
    also assumes that the files have the string 'gen' or 'mrna' in the filename,
    for genomic and mRNA sequence, respectively. Output is given as a series of
    csv values for various pertinent stats.""")
parser.add_argument('infiles', nargs='+', help='list of aligned infiles')
args = parser.parse_args()

def gulp(string, start, gulp_size):
    """get substrings of a string"""
    gulpstr = ''
    chars = string[start:start+gulp_size]
    for char in chars:
        gulpstr += char
    return gulpstr

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


for infile in args.infiles:
    name = infile.strip('.afa')
    with open(infile,'U') as f:
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
        if re.search('mrna',k) or re.search('mRNA',k) or re.search('mRNAs',k):
            mseq = seqdict.get(k)
        else:
            gseq = seqdict.get(k)

    smseq = sanitize(mseq)
    sgseq = sanitize(gseq)
    #print "smseq"
    #print smseq
    #print "sgseq"
    #print sgseq

    seq_pair = CodonPair(smseq, sgseq, name)  #sanitized sequences for direct comparison

    i = 0
    while not compare_seqs((gulp(mseq, i, 9)), (gulp(gseq, i, 9))):  #start of alignment
        if gseq[i] != '-':
            seq_pair.incr_all()
        if mseq[i] != '-':
            seq_pair.incr_mrna()
        i += 1
    j = 0
    while not compare_seqs((gulp(mseq[::-1], j, 9)), (gulp(gseq[::-1], j, 9))):  #end of alignment
        j += 1

    newmseq = mseq[i:(len(mseq)-j)]  #only compare regions that align
    newgseq = gseq[i:(len(gseq)-j)]
    #print "newmseq"
    #print newmseq
    #print "newgseq"
    #print newgseq

    for i, (rg, rm) in enumerate(zip(newgseq, newmseq)):  #compare matching regions
        #print i
        #print "codon position" + '' + str(seq_pair.codon_pos)
        if seq_pair.codon_pos != 3 or i < 2:  #only want codons once, only want full codons in both
            seq_pair.incr_all()
            seq_pair.incr_mrna()
        else:
            if rg == '-' and rm != '-':  #insertion in mRNA
                seq_pair.incr_mrna()
            elif rm == '-' and rg != '-':  #insertion in DNA
                seq_pair.incr_all()
            else:  #residue in both
                seq_pair.update_gcodons()
                seq_pair.update_mcodons()
                seq_pair.incr_all()
                seq_pair.incr_mrna()

    out = name + "_codons.csv"
    with open(out,'w') as o:
        o.write("amino acid,codon,genome usage,mRNA usage")
        o.write("\n")
        for k1, k2 in zip(seq_pair.gnuc_aa_dict, seq_pair.mnuc_aa_dict):  #loop through both dictionaries
            o.write(k1)
            for k3, k4 in zip(seq_pair.gnuc_aa_dict[k1].keys(), seq_pair.mnuc_aa_dict[k2].keys()):
                o.write(',' + k3 + ',' + str(seq_pair.gnuc_aa_dict[k1][k3]) + ',' + str(seq_pair.mnuc_aa_dict[k1][k4]) + "\n")
            o.write("\n")
