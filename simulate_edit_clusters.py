#!/usr/bin/env python

import os
import re
import sys
import argparse
import scipy.stats as st

from rediting.util import files, sequence, strings

parser = argparse.ArgumentParser(
    description = """Performs simulations to assess significance of edit clustering""",
    epilog = """This program is essentially identical to the sliding
    window program provided in the same package, except that it will
    require only aligned genomic/RNA sequence pairs, and will not
    output graphs of the calculated values. Instead, the program
    performs a number of simulations, in each case generating a sequence
    with the same length and number of edits as the actual genomic/
    RNA pair. The two distributions are compared to determine if they
    could have arisen from populations (sequences) with equal variance.""")
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
parser.add_argument('-x', '--simulations', help='number of simulations', default=100)
args = parser.parse_args()

window_size = float(args.window_size)
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

# Find the beginning and end of aligned region
i = 0
j = 0
try:
    # Compare genomic and RNA sequences to find local regions of good
    # similarity, this is taken as the start and end of aligned region
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

# Once we know the start and end, simply chop off everything else
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
num_total_edits = 0
compstr1 = ''
for i, (rg, rm) in enumerate(zip(new_gen_seq, new_rna_seq)):
    # Insertion in mrna
    if rg == '-' and rm != '-':
        compstr1 += str(0)
    # Insertion in genomic
    elif rm == '-' and rg != '-':
        compstr1 += str(0)
    # No edits
    elif rg == rm:
        compstr1 += str(0)
    # Edit
    elif rg != rm:
        # The value 1 here indicates an edit
        compstr1 += str(1)
        num_total_edits += 1
    else:
        pass

# Calculates percent edits for each window
edit_list = []
# Gets all full-length windows possible for length of aligned sequences
for start,end in sequence.get_indices(compstr1, window_size):
    try:
        edit_list.append(sequence.calc_percent(
            compstr1, start, end, window_size))
    # This error should not get thrown, but just in case
    except(ValueError,IndexError):
        print "Error detected while adding to edit_list"
        pass

length = len(new_rna_seq)
num_edits = num_total_edits
p_values = []
x = 0
while x < num_gens:
    seq1 = sequence.generate_start_sequence(length,bases)
    # Add mutations to simulate edits
    seq2 = sequence.mutate_sequence(seq1,num_edits,bases)
    sim_dist = sequence.calc_cluster_score(seq1,seq2)
    # Perform test, median-centred
    p_val = st.levene(edit_list,sim_dist)
    p_values.append(p_val[1])

    x += 1

sig_pvals = 0
for p_val in p_values:
    if p_val < 0.05:
        sig_pvals += 1
percent_sig = (float(sig_pvals)/len(p_values))

m_o.write("%s,%s,%s,%.2f" % (k,length,num_edits,percent_sig))
m_o.write("\n")

# Finally close the output file again
m_o.close()
