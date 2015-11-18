#!/usr/bin/env python

import re
import os
import sys
import argparse

from classes import classes,matrices
from util import files, rmath, sequence, strings

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

m_out = args.outfile
# appends if specified file already exists
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    m_o = open(m_out,'w')
    m_o.write("gene,GC before,GC after,nucleotide length,amino acid length,"
        "percent edits,percent edits in first two positions,"
        "percent amino acid edits,average edit score")
    if args.polyt:
        m_o.write(",fraction polyT before,fraction polyT after,fraction " +
            str(percent) + " percent polyT before,fraction " + str(percent)
            + " percent polyT after")
    m_o.write("\n" * 2)

b_out = name + "_basic_editing.csv"
if args.edits:
    e_out = name + "_editing_types.txt"
if args.codon:
    c_out = name + "_codon_preference.csv"

seqdict = {}
files.build_seqdict(args.infile,seqdict)

rna_string = str(args.RNA)
gen_string = str(args.genomic)
for k in seqdict.keys():
    if re.search(rna_string,k):
        rna_seq = seqdict.get(k)
    elif re.search(gen_string,k):
        gen_seq = seqdict.get(k)

# We directly compare aligned sequences, but class implementation uses
# unaligned sequences (i.e. no gap characters '-')
san_rna_seq = strings.sanitize(rna_seq)
san_gen_seq = strings.sanitize(gen_seq)
seq_pair = classes.SeqPair(san_rna_seq,san_gen_seq,name)

# Need to find beginning and end of aligned region
i = 0
j = 0
try:
    while not sequence.compare_seqs((strings.gulp(rna_seq, i, size)),
            (strings.gulp(gen_seq, i, size)), num_equal):
        if gen_seq[i] != '-':
            seq_pair.incr_all()
        if rna_seq[i] != '-':
            seq_pair.incr_mrna()
        i += 1
    while not sequence.compare_seqs((strings.gulp(rna_seq[::-1], j, size)),
            (strings.gulp(gen_seq[::-1], j, size)), num_equal):
        j += 1
# If we get an index error then we cannot find start and end of both sequences
except(IndexError):
    print "Could not discern aligned part of sequences for gene " + str(name)
    # Exit cleanly
    sys.exit(0)

new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]

edit_list = []
if args.edits:
    num_edited_res = 0

# Compare matching regions and look for unequal residues (edits)
for i, (rg, rm) in enumerate(zip(new_gen_seq, new_rna_seq)):
    # If there is an insertion in the RNA only
    if rg == '-' and rm != '-':
        seq_pair.incr_mrna()
    # If there is an insertion in the DNA only
    elif rm == '-' and rg != '-':
        seq_pair.incr_all()
    # Base present in both, but no editing
    elif rg == rm:
        if args.codon:
            # Saves us from updating same codon more than once
            if seq_pair.codon_pos != 3 or i < 2:
                pass
            else:
                seq_pair.update_gcodons()
                seq_pair.update_mcodons()
        seq_pair.incr_all()
        seq_pair.incr_mrna()
    # Base present in both, but there is an edit!
    elif rg != rm:
        pos = seq_pair.index_nuc() + 1
        cpos = seq_pair.index_position()
        gnuc = seq_pair.lookup_gnuc()
        mnuc = seq_pair.lookup_mnuc()
        gcod = seq_pair.lookup_gcodon()
        mcod = seq_pair.lookup_mcodon()
        gaa = seq_pair.lookup_gaa()
        maa = seq_pair.lookup_maa()
        scr = (matrices.Blosum62(gaa, maa).sub_score())

        if args.polyt:
            is_polyt = "N"
            # Only look at first seven bases at the start
            if i <= 4:
                polyt_test_seq = strings.gulp(new_gen_seq, 0, 7)
            # Only look at last seven bases at the end
            elif i >= len(new_gen_seq) - 4:
                polyt_test_seq = strings.gulp(new_gen_seq,
                        len(new_gen_seq)-7, 7)
            # In the middle take 3 bases on either side (seven total)
            else:
                polyt_test_seq = strings.gulp(new_gen_seq, i-3, 7)
            if sequence.polyT(polyt_test_seq):
                is_polyt = "Y"

            is_polyt_percent = "N"
            percent_polyt_seqs = []
            for y in range(10):
                try:
                    percent_polyt_seqs.append(
                            strings.gulp(new_gen_seq, i-y, 10))
                # Should only fail towards end of sequence
                except(IndexError):
                    pass
            if rmath.ispolyTpercent(percent_polyt_seqs, percent):
                is_polyt_percent = "Y"

        if not args.polyt:
            edit_list.append([pos,cpos,gnuc,mnuc,gcod,mcod,gaa,maa,scr])
        elif args.polyt:
            edit_list.append([pos,cpos,gnuc,mnuc,gcod,mcod,gaa,maa,scr,
                is_polyt,is_polyt_percent])

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

# Since polyT does not rely on presence/absence of edits per se, simplest
# thing to do is to run through each sequence again and count total occurrence
if args.polyt:
    num_gen_polyt = 0.0
    num_rna_polyt = 0.0
    polyt_indices = sequence.get_indices(new_gen_seq, 7)

    for start,end in polyt_indices:
        if sequence.polyT(new_gen_seq[start:end]):
            num_gen_polyt += 1.0
        if sequence.polyT(new_rna_seq[start:end]):
            num_rna_polyt += 1.0

    num_gen_percent_polyt = 0.0
    num_rna_percent_polyt = 0.0
    percent_polyt_indices = sequence.get_indices(new_gen_seq, 10)

    for start,end in percent_polyt_indices:
        if sequence.polyTpercent(new_gen_seq[start:end], percent):
            num_gen_percent_polyt += 1.0
        if sequence.polyTpercent(new_rna_seq[start:end], percent):
            num_rna_percent_polyt += 1.0

        # Calculate the fraction, since we take overlapping indices
        fraction_gen_polyt = (num_gen_polyt/len(polyt_indices)) * 100
        fraction_rna_polyt = (num_rna_polyt/len(polyt_indices)) * 100
        fraction_gen_percent_polyt = (num_gen_percent_polyt/
                len(percent_polyt_indices)) * 100
        fraction_rna_percent_polyt = (num_rna_percent_polyt/
                len(percent_polyt_indices)) * 100

    subscore = 0
    num_aaedits = 0
    num_fpos = 0
    with open(b_out,'w') as b_o:
        b_o.write("position,codon position,genome base,mRNAbase,genome codon,"
            "mRNA codon,genome amino acid,mRNA amino acid,substitution score")
        if args.polyt:
            b_o.write(",in a polyT tract,in a tract with " +
                    str(percent) + " percent T residues")
        b_o.write("\n" * 2)
        if args.polyt:
            for P, C, GN, MN, GC, MC, GA, MA, S, IP, IPP in edit_list:
                b_o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s"
                        % (P,C,GN,MN,GC,MC,GA,MA,S,IP,IPP) + "\n")
                subscore += int(S)
                if GA != MA:
                    num_aaedits += 1
                if C == 1 or C == 2:
                    num_fpos += 1
        else:
            for P, C, GN, MN, GC, MC, GA, MA, S in edit_list:
                b_o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s" %
                        (P,C,GN,MN,GC,MC,GA,MA,S) + "\n")
                subscore += int(S)
                if GA != MA:
                    num_aaedits += 1
                if C == 1 or C == 2:
                    num_fpos += 1

    gcb = sequence.calc_gc(new_gen_seq)
    gca = sequence.calc_gc(new_rna_seq)
    seqlength = len(new_rna_seq)
    aalength = seqlength/3
    numedits = float(len(edit_list))
    seqedits = (numedits/seqlength) * 100
    aaedits = (float(num_aaedits)/aalength) * 100
    try:
        editscore = subscore/numedits
    # If no residues are edited then this will throw an error
    except(ZeroDivisionError):
        editscore = 0
    # Again, same error will occur if numedits is zero
    try:
        fpos = (num_fpos/numedits) * 100
    except(ZeroDivisionError):
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
