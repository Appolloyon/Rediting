#!/usr/bin/env python

import os
import re
import sys
import argparse
import numpy as np

from rediting.util import files, sequence, strings, rmath

parser = argparse.ArgumentParser(
    description = """Calculates percent of edits occurring in regions
        with higher than average editing rate across the gene""",
    epilog = """This program is essentially identical to the sliding
    window program provided in the same package, except that it will
    require only aligned genomic/RNA sequence pairs, and will not
    output graphs of the calculated values. Instead, it will calculate
    the proportion of edits occurring in regions with higher than average
    editing, as well as a 'cluster score', which is related to the window
    value for the edit and the editing mean.""")
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
parser.add_argument('--simulation', action='store_true')
args = parser.parse_args()

window_size = float(args.window_size)
num_equal = int(args.numequal)
size = int(args.size)
name = args.name

if not args.simulation:
    # Create a "master" outfile to collate data from multiple files
    m_out = args.outfile
    # Appends if specified file already exists
    if os.path.isfile(m_out):
        m_o = open(m_out,'a')
    else:
        # The first time the file is opened, write header lines
        m_o = open(m_out,'w')
        m_o.write("name,total edits,edits above average,edits at or below average,"
            "percent above average edits,cluster score")
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

# Get the average of all edits over windows
# Note that this will likely differ slightly from
# a global calculation, i.e. value/len*100
try:
    edit_mean = rmath.calc_mean(edit_list)
    edit_std_dev = np.std(edit_list)
# If no edits occurred, then edit_list is empty
# and calc_mean will throw a ZeroDivError
except(ZeroDivisionError):
    print "Zero Div Error calculating edit_mean"
    edit_mean = 0.0
    edit_std_dev = 0.0

# Determine how many edits are above the average
edits_above_average = 0.0
total_edits = 0.0
cluster_score = 0.0
# Here we calculate edits based on the actual percent in each window
for edit in edit_list:
    if edit > edit_mean:
        total_edits += edit
        edits_above_average += edit
    else:
        total_edits += edit
    try:
        #cluster_score += (edit * (abs(edit - edit_mean)/edit_std_dev))
        cluster_score += (edit * (abs(edit - edit_mean)**2))
    except(ZeroDivisionError):
        cluster_score += 0.0


# Whatever is left becomes "other"
other_edits = total_edits - edits_above_average

# Determine percent of edits above average
try:
    percent_above_average_edits = (edits_above_average/total_edits) * 100
    percent_other_edits = (other_edits/total_edits) * 100
# Again, if no edits occured this will throw an error
except(ZeroDivisionError):
    print "Zero Div Error calculating above average edits"
    percent_above_average_edits = 0.0
    percent_other_edits = 0.0

# This seems like an odd way to calculate it, but we need to convert
# a percent over many windows back to a number
num_above_average = round((percent_above_average_edits/100) * num_total_edits)
num_other = round((percent_other_edits/100) * num_total_edits)
try:
    percent_diff = (num_above_average/num_total_edits) * 100
except(ZeroDivisionError):
    percent_diff = 0.0

if not args.simulation:
    m_o.write("%s,%s,%s,%s,%.2f,%.2f" % (name,num_total_edits,num_above_average,\
        num_other,percent_diff,cluster_score))
    m_o.write("\n")

if not args.simulation:
    # Finally close the output file again
    m_o.close()

if args.simulation:
    tmpfile = 'tempfile.csv'
    with open(tmpfile,'w') as o:
        o.write("%s,%s,%s,%f" % (name,len(new_rna_seq),num_total_edits,cluster_score))
        for edit in edit_list:
            o.write(",%f" % (edit))
        o.write("\n")
