#!/usr/bin/env python

import re
#import sys
import argparse

from classes import CodonPair
from matrices import Blosum62
from functions import gulp, compare_seqs, nonblank_lines, sanitize, calc_gc, build_seqdict

parser = argparse.ArgumentParser(
    description = """Calculates editing stats between genomic/RNA sequences""",
    epilog = """This program assumes that the genomic and RNA sequences for a
    given gene are provided in the same file in an FASTA format, such as that
    output by MAFFT or MUSCLE. It will go through and calculate information
    regarding the editing events, as identified by differences in the aligned
    sequences. These include nucleotide, codon, and amino acid changes, as well
    as summarizing changes in all of these as well as GC content. Depending on
    user input, one or more files will be created either in .txt or .csv format.
    In order to distinguish between genomic and RNA sequences, the user must
    specify a distinguishing string (word or list of characters) present in
    the FASTA header of each (e.g. 'RNA' or 'mRNA' for RNA sequences.""")
parser.add_argument('infiles', nargs='+', help='list of aligned infiles')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence header')
parser.add_argument('-g', '--genomic', help='unique string present in genomic sequence headers')
parser.add_argument('-n', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignment', default=7)
parser.add_argument('-s', '--size', help='number of residues to compare to determine start/end\
        of an alignment', default=9)
parser.add_argument('-b', '--basic', action='store_true', help='calculate only basic editing stats')
parser.add_argument('-e', '--edits', action='store_true', help='summarize editing types')
parser.add_argument('-c', '--codon', action='store_true', help='summarize codon usage difference')
args = parser.parse_args()

#if not args.basic and not args.edits and not args.codon:
#    print args.description
#    sys.exit("please specify one of at least basic, edits, or codon")

if args.basic:
    m_out = "master_editing_out.csv"
    m_o = open(m_out,'w')

for infile in args.infiles:
    name = infile.split('.')[0]
    gene = "pass"

    if args.basic:
        b_out = name + "_basic_editing.csv"
    if args.edits:
        e_out = name + "_editing_types.txt"
    if args.codon:
        c_out = name + "_codon_preference.csv"

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
    seq_pair = CodonPair(san_rna_seq,san_gen_seq,name)

    num_equal = int(args.numequal)
    size = int(args.size)
    i = 0
    j = 0
    while not compare_seqs((gulp(rna_seq, i, size)),\
            (gulp(gen_seq, i, size)), num_equal):
        if gen_seq[i] != '-':
            seq_pair.incr_all()
        if rna_seq[i] != '-':
            seq_pair.incr_mrna()
        i += 1
    while not compare_seqs((gulp(rna_seq[::-1], j, size)),\
            (gulp(gen_seq[::-1], j, size)), num_equal):
        j += 1

    new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
    new_gen_seq = gen_seq[i:(len(gen_seq)-j)]

    if args.basic:
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
            if args.basic:
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

    if args.basic:
        subscore = 0
        num_aaedits = 0
        num_fpos = 0
        with open(b_out,'w') as b_o:
            b_o.write("position,codon position,genome base,mRNA base,genome codon,\
               mRNA codon,genome amino acid,mRNA amino acid,substitution score")
            b_o.write("\n")
            for P, C, GN, MN, GC, MC, GA, MA, S in edit_list:
                b_o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s" % (P,C,GN,MN,GC,MC,GA,MA,S) + "\n")
                subscore += int(S)
                if GA != MA:
                    num_aaedits += 1
                    #print num_aaedits
                if C == 1 or C == 2:
                    num_fpos += 1

        gcb = calc_gc(new_gen_seq)
        gca = calc_gc(new_rna_seq)
        seqlength = len(new_rna_seq)
        aalength = seqlength/3
        numedits = float(len(edit_list))
        seqedits = (numedits/seqlength) * 100
        aaedits = (float(num_aaedits)/aalength) * 100
        editscore = subscore/numedits
        fpos = (num_fpos/numedits) * 100

        m_o.write("%s,%.2f,%.2f,%s,%s,%.2f,%.2f,%.2f,%.2f" % (gene,gcb,gca,seqlength,\
                aalength,seqedits,fpos,aaedits,editscore) + "\n")

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

"""
Program Outline:
    For each infile:
        Make all necessary outfiles
        Parse the infile and build the seqdict
        Identify which sequences are the mrna and genomic ones
        Make necessary objects
        Iterate through sequences and populate sequence objects as necessary
        Write all information to outfiles
    Done
"""
