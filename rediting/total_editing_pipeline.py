#!/usr/bin/env python

import re
import os
import sys
import argparse

from classes import SeqPair
from matrices import Blosum62
from functions import gulp, compare_seqs, sanitize, calc_gc, build_seqdict,\
        get_indices, polyT, polyTpercent, ispolyTpercent

parser = argparse.ArgumentParser(
    description = """Calculates editing stats between genomic/RNA sequences""",
    epilog = """This program assumes that the genomic and RNA sequences for a
    given gene are provided in the same file in an aligned FASTA format, such as that
    output by MAFFT or MUSCLE. It will go through and calculate information
    regarding the editing events, as identified by differences in the aligned
    sequences. These include nucleotide, codon, and amino acid changes, as well
    as summarizing changes in all of these as well as GC content. Depending on
    user input, one or more files will be created either in .txt or .csv format.
    In order to distinguish between genomic and RNA sequences, the user must
    specify a distinguishing string (word or list of characters) present in
    the FASTA header of each (e.g. 'RNA' or 'mRNA' for RNA sequences.""")
parser.add_argument('-in', '--infile', help='infile with aligned sequences')
parser.add_argument('-out', '--outfile', help='name for master outfile')
parser.add_argument('-n', '--name', help='name to append to output file')
parser.add_argument('-g', '--gene', help='gene name for output files')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence header')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequence headers')
parser.add_argument('-neq', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignment', default=7)
parser.add_argument('-s', '--size', help='number of residues to compare to determine start/end\
        of an alignment', default=9)
parser.add_argument('-e', '--edits', action='store_true', help='summarize editing types')
parser.add_argument('-c', '--codon', action='store_true', help='summarize codon usage difference')
parser.add_argument('-t', '--polyt', action='store_true', help='calculate polyT')
parser.add_argument('-p', '--percent', help='percent cut-off for polyT', default=70)
args = parser.parse_args()

num_equal = int(args.numequal)
size = int(args.size)
percent = int(args.percent)
name = args.name
gene = args.gene

m_out = args.outfile #"master_editing_out.csv"
# appends if specified file already exists
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    m_o = open(m_out,'w')
    m_o.write("gene,GC before,GC after,nucleotide length,amino acid length,percent edits,percent edits in first two positions,percent amino acid edits,average edit score")
    if args.polyt:
        m_o.write(",fraction polyT before,fraction polyT after,fraction " + str(percent) + " percent polyT before,fraction " + str(percent) + " percent polyT after")
    else:
        pass
    m_o.write("\n" * 2)

#for infile in args.infiles:
    #name = infile.split('.')[0]
    #gene = name.split('_')[1]

b_out = name + "_basic_editing.csv"
if args.edits:
    e_out = name + "_editing_types.txt"
if args.codon:
    c_out = name + "_codon_preference.csv"

seqdict = {}
build_seqdict(args.infile,seqdict)

rna_string = str(args.RNA)
gen_string = str(args.genomic)
for k in seqdict.keys():
    if re.search(rna_string,k):
        rna_seq = seqdict.get(k)
    elif re.search(gen_string,k):
        gen_seq = seqdict.get(k)

san_rna_seq = sanitize(rna_seq)
san_gen_seq = sanitize(gen_seq)
seq_pair = SeqPair(san_rna_seq,san_gen_seq,name)

i = 0
j = 0
try:
    while not compare_seqs((gulp(rna_seq, i, size)),
            (gulp(gen_seq, i, size)), num_equal):
        if gen_seq[i] != '-':
            seq_pair.incr_all()
        if rna_seq[i] != '-':
            seq_pair.incr_mrna()
        i += 1
    while not compare_seqs((gulp(rna_seq[::-1], j, size)),
            (gulp(gen_seq[::-1], j, size)), num_equal):
        j += 1
except(IndexError):
    print "Could not discern aligned part of sequences"
    sys.exit(0)

new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]

edit_list = []
if args.edits:
    num_edited_res = 0

for i, (rg, rm) in enumerate(zip(new_gen_seq, new_rna_seq)):  #compare matching regions
    if rg == '-' and rm != '-':  #insertion in mRNA
        seq_pair.incr_mrna()
    elif rm == '-' and rg != '-':  #insertion in DNA
        seq_pair.incr_all()
    elif rg == rm:  #residue in both, but no edits
        if args.codon:
            if seq_pair.codon_pos != 3 or i < 2:
                pass
            else:
                seq_pair.update_gcodons()
                seq_pair.update_mcodons()
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

        if args.polyt:
            is_polyt = "N"
            if i <= 4:
                polyt_test_seq = gulp(new_gen_seq, 0, 7)
            elif i >= len(new_gen_seq) - 4:
                polyt_test_seq = gulp(new_gen_seq, len(new_gen_seq)-7, 7)
            else:
                polyt_test_seq = gulp(new_gen_seq, i-3, 7)
            if polyT(polyt_test_seq):
                is_polyt = "Y"

            is_polyt_percent = "N"
            percent_polyt_seqs = []
            for y in range(10):
                try:
                    percent_polyt_seqs.append(gulp(new_gen_seq, i-y, 10))
                except:
                    pass
            if ispolyTpercent(percent_polyt_seqs, percent):
                is_polyt_percent = "Y"

        if not args.polyt: #args.basic and not args.polyt:
            edit_list.append([pos,cpos,gnuc,mnuc,gcod,mcod,gaa,maa,scr])
        elif args.polyt: #args.basic and args.polyt:
            edit_list.append([pos,cpos,gnuc,mnuc,gcod,mcod,gaa,maa,scr,is_polyt,is_polyt_percent])

        if args.codon:
            if seq_pair.codon_pos != 3 or i < 2:
                pass
            else:
                seq_pair.update_gcodons()
                seq_pair.update_mcodons()

        if args.edits:
            seq_pair.update_transdict()
            num_edited_res += 1

        seq_pair.incr_all()
        seq_pair.incr_mrna()

if args.polyt:
    num_gen_polyt = 0.0
    num_rna_polyt = 0.0
    polyt_indices = get_indices(new_gen_seq, 7)

    for start,end in polyt_indices:
        if polyT(new_gen_seq[start:end]):
            num_gen_polyt += 1.0
        if polyT(new_rna_seq[start:end]):
            num_rna_polyt += 1.0

    num_gen_percent_polyt = 0.0
    num_rna_percent_polyt = 0.0
    percent_polyt_indices = get_indices(new_gen_seq, 10)

    for start,end in percent_polyt_indices:
        if polyTpercent(new_gen_seq[start:end], percent):
            num_gen_percent_polyt += 1.0
        if polyTpercent(new_rna_seq[start:end], percent):
            num_rna_percent_polyt += 1.0

        fraction_gen_polyt = (num_gen_polyt/len(polyt_indices)) * 100
        fraction_rna_polyt = (num_rna_polyt/len(polyt_indices)) * 100
        fraction_gen_percent_polyt = (num_gen_percent_polyt/len(percent_polyt_indices)) * 100
        fraction_rna_percent_polyt = (num_rna_percent_polyt/len(percent_polyt_indices)) * 100

    subscore = 0
    num_aaedits = 0
    num_fpos = 0
    with open(b_out,'w') as b_o:
        b_o.write("position,codon position,genome base,mRNAbase,genome codon,mRNA codon,\
genome amino acid,mRNA amino acid,substitution score")
        if args.polyt:
            b_o.write(",in a polyT tract,in a tract with " + str(percent) + " percent T residues")
        b_o.write("\n" * 2)
        if args.polyt:
            for P, C, GN, MN, GC, MC, GA, MA, S, IP, IPP in edit_list:
                b_o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s" % (P,C,GN,MN,GC,MC,GA,MA,S,IP,IPP) + "\n")
                subscore += int(S)
                if GA != MA:
                    num_aaedits += 1
                if C == 1 or C == 2:
                    num_fpos += 1
        else:
            for P, C, GN, MN, GC, MC, GA, MA, S in edit_list:
                b_o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s" % (P,C,GN,MN,GC,MC,GA,MA,S) + "\n")
                subscore += int(S)
                if GA != MA:
                    num_aaedits += 1
                if C == 1 or C == 2:
                    num_fpos += 1

    gcb = calc_gc(new_gen_seq)
    gca = calc_gc(new_rna_seq)
    seqlength = len(new_rna_seq)
    aalength = seqlength/3
    numedits = float(len(edit_list))
    seqedits = (numedits/seqlength) * 100
    aaedits = (float(num_aaedits)/aalength) * 100
    try:
        editscore = subscore/numedits
    except:
        editscore = 0
    try:
        fpos = (num_fpos/numedits) * 100
    except:
        fpos = 0

    m_o.write("%s,%.2f,%.2f,%s,%s,%.2f,%.2f,%.2f,%.2f" % (gene,gcb,gca,seqlength,\
            aalength,seqedits,fpos,aaedits,editscore))
    if args.polyt:
        m_o.write(",%.2f,%.2f,%.2f,%.2f" % (fraction_gen_polyt,fraction_rna_polyt,\
            fraction_gen_percent_polyt,fraction_rna_percent_polyt))
    m_o.write("\n")

    if args.codon:
        with open(c_out,'w') as c_o:
            c_o.write("amino acid,codon,genome usage,mRNA usage")
            c_o.write("\n")
            # loop through both dictionaries
            for k1, k2 in zip(seq_pair.gnuc_aa_dict, seq_pair.mnuc_aa_dict):
                c_o.write(k1)
                for k3, k4 in zip(seq_pair.gnuc_aa_dict[k1].keys(),\
                        seq_pair.mnuc_aa_dict[k2].keys()):
                    c_o.write(',' + k3 + ',' + str(seq_pair.gnuc_aa_dict[k1][k3])\
                            + ',' + str(seq_pair.mnuc_aa_dict[k1][k4]) + "\n")
                c_o.write("\n")

    if args.edits:
        with open(e_out,'w') as e_o:
            e_o.write("Total number of unequal residues: {}\n".format(num_edited_res))
            e_o.write("Number of transitions: \n")
            e_o.write("A to T: {}\n".format(seq_pair.transition_dict.get('a_t')))
            e_o.write("A to G: {}\n".format(seq_pair.transition_dict.get('a_g')))
            e_o.write("A to C: {}\n".format(seq_pair.transition_dict.get('a_c')))
            e_o.write("T to A: {}\n".format(seq_pair.transition_dict.get('t_a')))
            e_o.write("T to G: {}\n".format(seq_pair.transition_dict.get('t_g')))
            e_o.write("T to C: {}\n".format(seq_pair.transition_dict.get('t_c')))
            e_o.write("G to A: {}\n".format(seq_pair.transition_dict.get('g_a')))
            e_o.write("G to T: {}\n".format(seq_pair.transition_dict.get('g_t')))
            e_o.write("G to C: {}\n".format(seq_pair.transition_dict.get('g_c')))
            e_o.write("C to A: {}\n".format(seq_pair.transition_dict.get('c_a')))
            e_o.write("C to T: {}\n".format(seq_pair.transition_dict.get('c_t')))
            e_o.write("C to G: {}\n".format(seq_pair.transition_dict.get('c_g')))

m_o.close()
