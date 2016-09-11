#!/usr/bin/env python

import os
import re
import sys
import argparse

from rediting.util import files, sequence, strings

parser = argparse.ArgumentParser(
    description = """Calculates amino acids before and after editing to a reference""",
    epilog = """Unlike the other programs in this suite that take aligned nucleotide
        sequences as input, this program is designed to work with aligned amino acid
        sequences instead. It will simply go through each sequence and determine
        whether the amino acid changes observed between genomic and RNA sequences
        make the resulting sequence more or less similar to the reference sequence.""")
parser.add_argument('-in', '--infile', help='infile with aligned sequences')
parser.add_argument('-out', '--outfile', help='name for master outfile')
parser.add_argument('-n', '--name', help='name to append to outfile file')
parser.add_argument('-g', '--gene', help='gene name for output files')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence headers')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequence headers')
parser.add_argument('-neq', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignment', default=5)
parser.add_argument('-s', '--size', help='number of residues to compare to determine start/end\
        of an alignment', default=10)
args = parser.parse_args()

num_equal = int(args.numequal)
size = int(args.size)
name = args.name
gene = args.gene

# Create a "master" outfile to collate data from multiple files
m_out = args.outfile
# Appends if specified file already exists
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    # The first time the file is opened, write header lines
    m_o = open(m_out,'w')
    m_o.write("gene,num AA changes,num identical before,num identical after,"
        "num similar before,num similar after,avgerage edit score diff")
    m_o.write("\n" * 2)

# Each independent file is also used to create its own file
b_out = name + "_aminoacid_changes.csv"

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

# Need to find beginning and end of aligned region
i = 0
j = 0
# Need to keep track of gen index
gen_index = 0
try:
    # Compare genomic and RNA sequences to find local regions of good
    # similarity, this is taken as the start and end of aligned region
    while not sequence.compare_seqs((strings.gulp(rna_seq, i, size)),
            (strings.gulp(gen_seq, i, size)), num_equal):
        if gen_seq[i] != '-':
            # If there is actually an amino acid in the genomic sequence
            # we need to move the index ahead
            gen_index += 1
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
new_ref_seq = ref_seq[i:(len(ref_seq)-j)]

edit_list = []
sequence.compare_aa_seqs(gen_index,new_gen_seq,new_rna_seq,
        new_ref_seq,edit_list)

num_ident_before = 0
num_ident_after = 0
num_similar_before = 0
num_similar_after = 0
diff_score = 0
with open(b_out,'w') as b_o:
    b_o.write("AA position,reference AA,genomic AA,transcript AA,transcript AA"
        "identical before,identical after,genomic AA similar,transcript AA similar,"
        "edit score difference")
    b_o.write("\n" * 2)
    for P, RA, GA, MA, IB, IA, SB, SA, D in edit_list:
        b_o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s" %
            (P,RA,GA,MA,IB,IA,SB,SA,D) + "\n")
        # Sum up all score differences to calculate and average later
        diff_score += int(D)
        # Keep a running tally on global values
        if IB == 'Y':
            num_ident_before += 1
        if IA == 'Y':
            num_ident_after += 1
        if SB == 'Y':
            num_similar_before += 1
        if SA == 'Y':
            num_similar_after += 1

# Each edit is one entry in the list of lists
numedits = float(len(edit_list))
try:
    avg_diff = diff_score/numedits
# If there are no edits, then this will throw an error
except(ZeroDivisionError):
    avg_diff = 0

m_o.write("%s,%s,%s,%s,%s,%s,%.2f" % (gene,numedits,num_ident_before,\
        num_ident_after,num_similar_before,num_similar_after,avg_diff))
m_o.write("\n")

m_o.close()
