#!/usr/bin/env python

import os
import re
import sys
import argparse

from util import files, sequence, strings, rmath

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
parser.add_argument('-w', '--window_size', help='size of sliding window', default=60)
args = parser.parse_args()

window_size = float(args.window_size)
num_equal = int(args.numequal)
size = int(args.size)
name = args.name

m_out = args.outfile
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    m_o = open(m_out,'w')
    m_o.write("name,average edit rate,percent above average edits")
    m_o.write("\n" * 2)

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

# Find the beginning and end of aligned region
i = 0
j = 0
try:
    while not sequence.compare_seqs((strings.gulp(rna_seq, i, size)),
            (strings.gulp(gen_seq, i, size)), num_equal):
        i += 1
    while not sequence.compare_seqs((strings.gulp(rna_seq[::-1], j, size)),
            (strings.gulp(gen_seq[::-1], j, size)), num_equal):
        j += 1
# If we get an index error then we cannot find start and end of both sequences
except(IndexError):
    print "Could not discern aligned part of sequences"
    # Exit cleanly
    sys.exit(0)

new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]

# Cannot do sliding window calculations if sequences are smaller than chosen
# window size; causes numerous errors in later parts of the program
if (len(new_rna_seq) < window_size or
    len(new_gen_seq) < window_size):
    print "One or more sequences are shorter than chosen window size."
    print "Please choose a smaller window size or check sequences."
    # Exit cleanly
    sys.exit(0)

# Calculate % edits between genomic and RNA sequences
compstr1 = ''
for i, (rg, rm) in enumerate(zip(new_gen_seq, new_rna_seq)):
    if rg == '-' and rm != '-':
        compstr1 += str(0)
    elif rm == '-' and rg != '-':
        compstr1 += str(0)
    elif rg == rm:
        compstr1 += str(0)
    elif rg != rm:
        compstr1 += str(1)
    else:
        pass

# Calculates percent edits for each window
edit_list = []
for start,end in sequence.get_indices(compstr1, window_size):
    try:
        edit_list.append(sequence.calc_percent(
            compstr1, start, end, window_size))
    # This error should not get thrown, but just in case
    except(ValueError,IndexError):
        print "Error detected while adding to edit_list"
        pass

# Get the average of all edits over windows
try:
    edit_mean = rmath.calc_mean(edit_list)
# If no edits occurred, then edit_list is empty
# and calc_mean will throw a ZeroDivError
except(ZeroDivisionError):
    print "Zero Div Error calculating edit_mean"
    edit_mean = 0.0

# Determine how many edits are above the average
edits_above_average = 0.0
total_edits = 0.0
for edit in edit_list:
    if edit > edit_mean:
        total_edits += edit
        edits_above_average += edit
    else:
        total_edits += edit

# Determine percent of edits above average
try:
    percent_above_average_edits = (edits_above_average/total_edits) * 100
# Again, if no edits occured this will throw an error
except(ZeroDivisionError):
    print "Zero Div Error calculating above average edits"
    percent_above_average_edits = 0.0

m_o.write("%s,%.2f,%.2f" % (name,edit_mean,percent_above_average_edits))
m_o.write("\n")

# Finally close the output file again
m_o.close()
