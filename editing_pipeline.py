#!/usr/bin/env python

import re
import sys
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
    while not compare_seqs((gulp(rna_seq, i, size)), (gulp(gen_seq, i, size)), num_equal):
        if gen_seq[i] != '-':
            seq_pair.incr_all()
        if rna_seq[i] != '-':
            seq_pair.incr_mrna()
        i += 1
    while not compare_seqs((gulp(rna_seq[::-1], j, size)), (gulp(gen_seq[::-1], j, size)), num_equal):
        j += 1

    new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
    new_gen_seq = gen_seq[i:(len(gen_seq)-j)]





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
