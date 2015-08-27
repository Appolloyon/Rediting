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
from classes import SeqPair

parser = argparse.ArgumentParser(
    description = """Calculates RNA editing in genomic poly-T tracks""",
    epilog = """This program assumes that the genomic and RNA sequence for a
    given gene are provided in the same file in aligned FASTA format (.afa). It
    also assumes that the files have the string 'gen' or 'mrna' in the seqname,
    for genomic and mRNA sequence, respectively. Output is given as a series of
    csv values for various pertinent stats.""")
parser.add_argument('-p', '--percent', help='percent cutoff for polyT')
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

def ispolyTpercent(plist):
    """check list elements for at least one polyT stretch"""
    for e in plist:
        if polyTpercent(e, percent):
            return True
        else:
            pass
    return False


for infile in args.infiles:
    name = infile.strip('.afa')
    gene = name.split('_')[1]
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

    seq_pair = SeqPair(smseq, sgseq, name)  #sanitized sequences for direct comparison
    seq_pair2 = SeqPair(smseq, sgseq, name)

    i = 0
    while not compare_seqs((gulp(mseq, i, 9)), (gulp(gseq, i, 9))):  #start of alignment
        if gseq[i] != '-':
            seq_pair.incr_all()
            seq_pair2.incr_all()
        if mseq[i] != '-':
            seq_pair.incr_mrna()
            seq_pair2.incr_mrna()
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

    edit_list = []

    for i, (rg, rm) in enumerate(zip(newgseq, newmseq)):  #compare matching regions
        if i <= 4:
            testseq = gulp(newgseq, 0, 7)
        elif i >= len(newgseq) - 4:
            testseq = gulp(newgseq, len(newgseq)-7, 7)
        else:
            testseq = gulp(newgseq, i-3, 7)
        if rg == '-' and rm != '-':  #insertion in mRNA
            seq_pair.incr_mrna()
        elif rm == '-' and rg != '-':  #insertion in DNA
            seq_pair.incr_all()
        elif rg == rm:  #residue in both, but no edits
            seq_pair.incr_all()
            seq_pair.incr_mrna()
        elif rg != rm:  #residue in both, but edited
            if polyT(testseq):
                pos = seq_pair.index_nuc() + 1
                cpos = seq_pair.index_position()
                gnuc = seq_pair.lookup_gnuc()
                mnuc = seq_pair.lookup_mnuc()
                gcod = seq_pair.lookup_gcodon()
                mcod = seq_pair.lookup_mcodon()
                gaa = seq_pair.lookup_gaa()
                maa = seq_pair.lookup_maa()
                edit_list.append([pos,cpos,gnuc,mnuc,gcod,mcod,gaa,maa])
                seq_pair.incr_all()
                seq_pair.incr_mrna()
            else:
                seq_pair.incr_all()
                seq_pair.incr_mrna()

    edit_dict2 = {}
    percent = int(args.percent)
    for i, (rg, rm) in enumerate(zip(newgseq, newmseq)):
        if rg == '-' and rm != '-':  #insertion in mRNA
            seq_pair2.incr_mrna()
        elif rm == '-' and rg != '-':  #insertion in DNA
            seq_pair2.incr_all()
        elif rg == rm:  #residue in both, but no edits
            seq_pair2.incr_all()
            seq_pair2.incr_mrna()
        elif rg != rm:  #residue in both, but edited
            testseq2 = []
            for y in range(10):
                try:
                    testseq2.append(gulp(newgseq, y-i, 10))
                except:
                    pass
            if ispolyTpercent(testseq2):
                #print i
                #print testseq2
                pos = seq_pair2.index_nuc() + 1
                cpos = seq_pair2.index_position()
                gnuc = seq_pair2.lookup_gnuc()
                mnuc = seq_pair2.lookup_mnuc()
                gcod = seq_pair2.lookup_gcodon()
                mcod = seq_pair2.lookup_mcodon()
                gaa = seq_pair2.lookup_gaa()
                maa = seq_pair2.lookup_maa()

                if pos not in edit_dict2.keys():
                    edit_dict2[pos] = []
                    edit_dict2[pos].extend([cpos,gnuc,mnuc,gcod,mcod,gaa,maa])

                seq_pair2.incr_all()
                seq_pair2.incr_mrna()
            else:
                seq_pair2.incr_all()
                seq_pair2.incr_mrna()

    out1 = name + "_polyt_out.csv"
    with open(out1,'w') as o1:
        o1.write("position,codon position,genome base,mRNA base,genome codon,\
mRNA codon,genome amino acid,mRNA amino acid")
        o1.write("\n")
        for P, C, GN, MN, GC, MC, GA, MA in edit_list:
            o1.write("%s,%s,%s,%s,%s,%s,%s,%s" % (P,C,GN,MN,GC,MC,GA,MA) + "\n")

    #print edit_list
    #print edit_dict2

    out2 = name + "_polyt_" + str(percent) + "%_out.csv"
    with open(out2,'w') as o2:
        o2.write("position,codon position,genome base,mRNA base,genome codon,\
mRNA codon,genome amino acid,mRNA amino acid")
        o2.write("\n")
        for pos in sorted(edit_dict2.keys()):
            o2.write("%s," % (pos))
            for elem in edit_dict2[pos]:
                #print elem
                o2.write("%s," % (elem))
            o2.write("\n")