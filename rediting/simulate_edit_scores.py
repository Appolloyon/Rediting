#!/usr/bin/env python

import os
import re
import sys
import argparse
import subprocess

from classes import classes, matrices
from util import files, sequence, strings

parser = argparse.ArgumentParser(
    description = """Calculates percent of edits occurring in regions
        with higher than average editing rate across the gene""",
    epilog = """This program is essentially identical to the sliding
    window program provided in the same package, except that it will
    require only aligned genomic/RNA sequence pairs, and will not
    output graphs of the calculated values.""")
parser.add_argument('-in', '--infile', help='infile with aligned sequences')
parser.add_argument('-out', '--outfile', help='name for master outfile')
parser.add_argument('-n', '--name', help='name to append to output file')
parser.add_argument('-r', '---RNA', help='unique string present in RNA sequence header')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequence header')
parser.add_argument('-neq', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignments', default=7)
parser.add_argument('-s', '--size', help='number of residues to compare to determine\
        start/end of an alignment', default=9)
parser.add_argument('-x', '--simulations', help='number of simulations', default=100)
args = parser.parse_args()

num_equal = int(args.numequal)
num_gens = int(args.simulations)
size = int(args.size)
name = args.name

# Global variable for possible bases in simulation
bases = 'AGTC'

# Create a "master" outfile to collate data from multiple files
m_out = args.outfile
# Appends if specified file already exists
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    # The first time the file is opened, write header lines
    m_o = open(m_out,'w')
    m_o.write("name,length,number edits,frequency of significant edits")
    m_o.write("\n" * 2)

# Load sequence data into a data structure for internal use
seqdict = {}
files.build_seqdict(args.infile,seqdict)

rna_string = str(args.RNA)
gen_string = str(args.genomic)
# Sequences must be in upper-case
for k in seqdict.keys():
    if re.search(rna_string,k):
        rna_seq = seqdict.get(k).upper()
    elif re.search(gen_string,k):
        gen_seq = seqdict.get(k).upper()
    else:
        ref_seq = seqdict.get(k).upper()

san_rna_seq = strings.sanitize(rna_seq)
san_gen_seq = strings.sanitize(gen_seq)
seq_pair = classes.SeqPair(san_rna_seq,san_gen_seq,name)

# Find the beginning and end of aligned region
i = 0
j = 0
try:
    # Compare genomic and RNA sequences to find local regions of good
    # similarity, this is taken as the start and end of aligned region
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
    print "Could not discern aligned part of sequences"
    # Exit cleanly
    sys.exit(0)

# Once we know the start and end, simply chop off everything else
new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]

num_edits,codon_positions,cumul_weights = sequence.get_positional_information(
        new_gen_seq,new_rna_seq,seq_pair.index_position())

with open("tempfile1.fa",'w') as o1:
    o1.write(">gen_seq" + "\n" + sequence.translate(gen_seq,1) + "\n")
    o1.write(">rna_seq" + "\n" + sequence.translate(rna_seq,1) + "\n")
    o1.write(">ref_seq" + "\n" + sequence.translate(ref_seq,1) + "\n")

subprocess.call(["muscle3.8.31_i86darwin64", "-in", "tempfile1.fa",
    "-out","tempfile1.afa"])


simulation = classes.Simulation()
gen_list = [b for b in gen_seq]
for i in range(num_gens):
    new_seq = "".join(sequence.weighted_mutation_new(gen_list,num_edits,
        codon_positions,cumul_weights,simulation))
    rna_start = rna_seq[:i]
    rna_end = rna_seq[(len(rna_seq)-j):]
    edited_rna_seq = rna_start + new_seq + rna_end

    with open("tempfile2.fa",'w') as o:
        o.write(">gen_seq" + "\n" + gen_seq)
        o.write(">rna_seq" + "\n" + edited_rna_seq)
        o.write(">ref_seq" + "\n" + ref_seq)
