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
from functions import gulp, compare_seqs, sanitize, build_seqdict, get_indices, polyT, ispolyTpercent

parser = argparse.ArgumentParser(
    description = """Calculates RNA editing in genomic poly-T tracks""",
    epilog = """This program assumes that the genomic and RNA sequence for a
    given gene are provided in the same file in aligned FASTA format (.afa). It
    also assumes that the files have the string 'gen' or 'mrna' in the seqname,
    for genomic and mRNA sequence, respectively. Output is given as a series of
    csv values for various pertinent stats.""")
parser.add_argument('-p', '--percent', help='percent cutoff for polyT')
parser.add_argument('infiles', nargs='+', help='list of aligned infiles')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence header')
parser.add_argument('-g', '--genomic', help='unique string present in genomic sequence header')
parser.add_argument('-n', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignment', default=7)
parser.add_argument('-s', '--size', help='number of residues to compare to determine start/end\
        of an alignment', default=9)
args = parser.parse_args()

for infile in args.infiles:
    name = infile.split('.')[0]
    gene = name.split('_')[1]

    seqdict = {}
    build_seqdict(infile,seqdict)

    rna_string = str(args.RNA)
    gen_string = str(args.genomic)
    for k in seqdict.keys():
        if re.search(rna_string,k):
            rna_seq = seqdict.get(k)
        elif re.search(gen_string,k):
            gen_seq = seqdict.get(k)

    san_rna_seq = sanitize(rna_seq)
    san_gen_seq = sanitize(gen_seq)

    seq_pair = SeqPair(san_rna_seq, san_gen_seq, name)  #sanitized sequences for direct comparison
    seq_pair2 = SeqPair(san_rna_seq, san_gen_seq, name)

    num_equal = int(args.numequal)
    size = int(args.size)
    i = 0
    j = 0
    while not compare_seqs((gulp(rna_seq, i, size)),
            (gulp(gen_seq, i, size)), num_equal):  #start of alignment
        if gen_seq[i] != '-':
            seq_pair.incr_all()
            seq_pair2.incr_all()
        if rna_seq[i] != '-':
            seq_pair.incr_mrna()
            seq_pair2.incr_mrna()
        i += 1
    while not compare_seqs((gulp(rna_seq[::-1], j, size)),
            (gulp(gen_seq[::-1], j, size)), num_equal):  #end of alignment
        j += 1

    newmseq = rna_seq[i:(len(rna_seq)-j)]  #only compare regions that align
    newgseq = gen_seq[i:(len(gen_seq)-j)]

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
                    testseq2.append(gulp(newgseq, i-y, 10))
                except:
                    pass
            #print i
            #print testseq2
            if ispolyTpercent(testseq2, int(args.percent)):
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

    num_gen_polyt = 0.0
    num_rna_polyt = 0.0
    num_gen_percent_polyt = 0.0
    num_rna_percent_polyt = 0.0

    for start,end in get_indices(newgseq, 7):
        gen_polyt_str = newgseq[start:end]
        rna_polyt_str = newmseq[start:end]
        #print "genomic sequence: " + gen_polyt_str
        #print "RNA sequence:     " + rna_polyt_str
        if polyT(gen_polyt_str):
            print "genomic sequence " + gen_polyt_str + " is a genomic polyT tract"
            num_gen_polyt += 1.0
        if polyT(rna_polyt_str):
            print "RNA sequence " + rna_polyt_str + " is a RNA polyT tract"
            num_rna_polyt += 1.0
    #for start,end in get_indices(newmseq, 7):
        #rna_polyt_str = newmseq[start:end]
        #print "rna sequence: " + rna_polyt_str
        #if polyT(rna_polyt_str):
            #print "is a polyT tract"
            #num_rna_polyt += 1.0
    percent_gen_polyt = (num_gen_polyt/len(get_indices(newgseq, 7))) * 100
    percent_rna_polyt = (num_rna_polyt/len(get_indices(newmseq, 7))) * 100

    print percent_gen_polyt
    print percent_rna_polyt

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
