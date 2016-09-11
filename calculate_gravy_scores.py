#!/usr/bin/env python

import re
import os
import sys
import argparse

from Bio.SeqUtils import ProtParam as PP
from rediting.util import files, sequence, strings

parser = argparse.ArgumentParser(
    description = """Calculates gravy scores before and after editing""",
    epilog = """This program takes aligned amino acid sequences for a
    genomic sequence and its edited RNA counterpart. It discerns the best
    aligned region and then calculates GRAVY hydrophobicity scores over
    this region. Positions corresponding to STOP signals and other non-
    standard amino acids are removed prior to calculating GRAVY scores.""")
parser.add_argument('-in', '--infile', help='infile with aligned sequences')
parser.add_argument('-out', '--outfile', help='name for master outfile')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence headers')
parser.add_argument('-g', '--gene', help='gene name for output file')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequence headers')
parser.add_argument('-neq', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignment', default=5)
parser.add_argument('-s', '--size', help='number of residues to compare to determine start/end\
        of an alignment', default=10)
args = parser.parse_args()

num_equal = int(args.numequal)
size = int(args.size)
gene = args.gene

# Define the standard amino acid alphabet
legal_aa = "GPAVLIMCFYWHKRQNEDST"

# Create a "master" outfile to collate data from multiple files
m_out = args.outfile
# Appends if specified file already exists
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    # The first time the file is opened, write header lines
    m_o = open(m_out,'w')
    m_o.write("gene,gen mol weight,gen gravy,rna mol weight,"
        "rna gravy")
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

# Need to find beginning and end of aligned region
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

# Need to remove positions in both sequences corresponding to
# internal gaps, STOP codons, or other non-legal amino acids
gen_indices = []
rna_indices = []
for (i,char) in enumerate(new_gen_seq):
    if char not in legal_aa:
        gen_indices.append(i)
for (i,char) in enumerate(new_rna_seq):
    if char not in legal_aa:
        rna_indices.append(i)
indices = list(set(gen_indices) | set(rna_indices))

# Now make the relevant BioPython objects
new_gen_seq = PP.ProteinAnalysis(sequence.trim_sequence(new_gen_seq,indices))
new_rna_seq = PP.ProteinAnalysis(sequence.trim_sequence(new_rna_seq,indices))

# Calculate parameters of interest for genomic seqs
gen_gravy = new_gen_seq.gravy()
gen_mol_weight = new_gen_seq.molecular_weight()
# Repeat for rna seqs
rna_gravy = new_rna_seq.gravy()
rna_mol_weight = new_rna_seq.molecular_weight()

m_o.write("%s,%.2f,%.2f,%.2f,%.2f" % (gene,gen_mol_weight,gen_gravy,\
        rna_mol_weight,rna_gravy))
m_o.write("\n")
