#!/usr/bin/env python

import re
import sys
import argparse

from rediting.util import files, sequence, strings, rmath

parser = argparse.ArgumentParser(
    description = """Calculates entropy scores compared to multiple references""",
    epilog = """Unlike the other programs in this suite that take aligned nucleotide
        sequences as input, this program is designed to work with aligned amino acid
        sequences instead. At each position in the alignment, as long as more than 50%
        of the reference sequences contain residues (not gaps), this script will calculate
        the positional entropy at the position in regards to all reference sequences.
        This is performed for both edited and unedited residues, and in the output file
        this is tracked in a binary way, i.e. 0 is not edited and 1 is edited.""")
parser.add_argument('-in', '--infile', help='infile with aligned sequences')
parser.add_argument('-n', '--name', help='name to append to outfile file')
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

# Each independent file is used to create its own file
b_out = name + "_entropy.csv"

# Load sequence data into a data structure for internal use
seqdict = {}
files.build_seqdict(args.infile,seqdict)

rna_string = str(args.RNA)
gen_string = str(args.genomic)
ref_keys = []
# Sequences must be in upper-case
for k in seqdict.keys():
    if re.search(rna_string,k):
        rna_seq = seqdict.get(k).upper()
    elif re.search(gen_string,k):
        gen_seq = seqdict.get(k).upper()
    else:
        # We want a list of all other headers
        ref_keys.append(k)

# Convert to counts by position in alignment
pos_dict = {}
for k,v in seqdict.items():
    if k in ref_keys: # We don't want to include the genomic/transcript sequences
        for i,aa in enumerate(v):
            if pos_dict.has_key((i)):
                pos_dict[(i)].append(aa)
            else:
                pos_dict[(i)] = []
                pos_dict[(i)].append(aa)

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

pos_list = []
# Calculates amino acid edits and entropy simultaeneously
for n, (rg,rm) in enumerate(zip(new_gen_seq, new_rna_seq)):
    # Insert in mrna sequence
    if rg == '-' and rm != '-':
        pass
    # Insert in genomic sequence
    elif rg != '-' and rm == '-':
        # Move the genomic index along
        gen_index += 1
    # Gaps in both sequences, no edits possible
    elif rg == '-' and rm == '-':
        pass
    # Residue in both, calculate entropy and move the genomic index along
    elif rg == rm:
        if sequence.check_gap_percent(pos_dict[i]):
            pos_entropy = rmath.calculate_entropy(pos_dict[i])
            pos_list.append([(i+1),(gen_index+1),0,pos_entropy,"0"])
        gen_index += 1
    # Edit detected, calculate entropy and move the genomic index along
    elif rg != rm:
        if sequence.check_gap_percent(pos_dict[i]):
            pos_entropy = rmath.calculate_entropy(pos_dict[i])
            avg_score = rmath.calculate_average_score_diff(rg,rm,pos_dict[i])
            pos_list.append([(i+1),(gen_index+1),1,pos_entropy,avg_score])
        gen_index += 1
    # No matter what we move through the alignment
    i += 1

with open(b_out,'w') as o:
    o.write("alignment position,genomic position,edited or not,positional entropy,average score difference")
    o.write("\n")
    for AP,GP,E,PE,S in pos_list:
        o.write("%s,%s,%d,%.2f," % (AP,GP,E,PE))
        try:
            o.write("%.2f" % (S))
        except TypeError,ValueError:
            o.write("%s" % (S))
        o.write("\n")
