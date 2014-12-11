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
from classes import SeqPair
from matrices import Blosum62

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


with open("m_out.csv",'w') as o2:
    o2.write("gene,GC before,GC after,nuc length,aa length,percent seq edits,\
percent edits in first two positions,percent aa edits,average edit score" + "\n")
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

        edit_list = []
        for i, (rg, rm) in enumerate(zip(newgseq, newmseq)):  #compare matching regions
            if rg == '-' and rm != '-':  #insertion in mRNA
                seq_pair.incr_mrna()
            elif rm == '-' and rg != '-':  #insertion in DNA
                seq_pair.incr_all()
            elif rg == rm:  #residue in both, but no edits
                seq_pair.incr_all()
                seq_pair.incr_mrna()
            elif rg != rm:  #residue in both, but edited
                pos = seq_pair.index_nuc() + 1
                cpos = seq_pair.index_position()
                gnuc = seq_pair.lookup_gnuc()
                mnuc = seq_pair.lookup_mnuc()
                gcod = seq_pair.lookup_gcodon()
                mcod = seq_pair.lookup_mcodon()
                gaa = seq_pair.lookup_gaa()
                maa = seq_pair.lookup_maa()
                scr = (Blosum62(gaa, maa).sub_score())
                edit_list.append([pos,cpos,gnuc,mnuc,gcod,mcod,gaa,maa,scr])
                seq_pair.incr_all()
                seq_pair.incr_mrna()

        out1 = name + "_out.csv"
        subscore = 0
        num_aaedits = 0
        num_fpos = 0
        with open(out1,'w') as o1:
            o1.write("position,codon position,genome base,mRNA base,genome codon,\
mRNA codon,genome amino acid,mRNA amino acid,substitution score")
            o1.write("\n")
            for P, C, GN, MN, GC, MC, GA, MA, S in edit_list:
                o1.write("%s,%s,%s,%s,%s,%s,%s,%s,%s" % (P,C,GN,MN,GC,MC,GA,MA,S) + "\n")
                subscore += int(S)
                if GA != MA:
                    num_aaedits += 1
                    #print num_aaedits
                if C == 1 or C == 2:
                    num_fpos += 1

        gcb = calc_gc(newgseq)
        gca = calc_gc(newmseq)
        seqlength = len(newmseq)
        aalength = seqlength/3
        numedits = float(len(edit_list))
        seqedits = (numedits/seqlength) * 100
        aaedits = (float(num_aaedits)/aalength) * 100
        editscore = subscore/numedits
        fpos = (num_fpos/numedits) * 100

        o2.write("%s,%.2f,%.2f,%s,%s,%.2f,%.2f,%.2f,%.2f" % (gene,gcb,gca,seqlength,\
                aalength,seqedits,fpos,aaedits,editscore) + "\n")
