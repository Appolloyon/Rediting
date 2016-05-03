#!/usr/bin/env python

import os
import re
import sys
import argparse
import subprocess
import scipy.stats as st

from classes import classes
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
parser.add_argument('-x', '--simulations', help='number of simulations', default=100)
args = parser.parse_args()

num_equal = int(args.numequal)
num_gens = int(args.simulations)
size = int(args.size)
name = args.name

# Create a "master" outfile to collate data from multiple files
m_out = args.outfile
# Appends if specified file already exists
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    # The first time the file is opened, write header lines
    m_o = open(m_out,'w')
    m_o.write("name,number edits,average edit score,average sim score,frequency of significant edits")
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
rna_start = rna_seq[:i]
rna_end = rna_seq[(len(rna_seq)-j):]
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]

num_edits,codon_positions,cumul_weights = sequence.get_positional_information(
        new_gen_seq,new_rna_seq,seq_pair.index_position())

with open("tempfile.fa",'w') as o1:
    o1.write(">gen_seq" + "\n" + sequence.translate(gen_seq,1) + "\n")
    o1.write(">mrna_seq" + "\n" + sequence.translate(rna_seq,1) + "\n")
    o1.write(">ref_seq" + "\n" + sequence.translate(ref_seq,1) + "\n")

subprocess.call(["muscle3.8.31_i86darwin64", "-in", "tempfile.fa",
    "-out", "tempfile.afa", "-quiet"])

aa_seqdict = {}
files.build_seqdict("tempfile.afa",aa_seqdict)
for k in aa_seqdict.keys():
    if re.search(rna_string,k):
        aa_rna_seq = aa_seqdict.get(k).upper()
    elif re.search(gen_string,k):
        aa_gen_seq = aa_seqdict.get(k).upper()
    else:
        aa_ref_seq = aa_seqdict.get(k).upper()

edit_list = []
sequence.compare_aa_seqs(0,aa_gen_seq,aa_rna_seq,aa_ref_seq,edit_list)
avg_score = rmath.calculate_score_diff(edit_list)
scr_diffs = []
for P,RA,GA,MA,IB,IA,SB,SA,D in edit_list:
    scr_diffs.append(int(D))

p_values = []
simulation = classes.Simulation()
total_sim_score = 0
gen_list = [b for b in new_gen_seq]
for i in range(num_gens):
    new_seq = "".join(sequence.weighted_mutation_new(gen_list,num_edits,
        codon_positions,cumul_weights,simulation))
    edited_rna_seq = rna_start + new_seq + rna_end

    with open("sim_tempfile.fa",'w') as o2:
        o2.write(">gen_seq" + "\n" + sequence.translate(gen_seq,1) + "\n")
        o2.write(">mrna_seq" + "\n" + sequence.translate(edited_rna_seq,1) + "\n")
        o2.write(">ref_seq" + "\n" + sequence.translate(ref_seq,1) + "\n")

    subprocess.call(["muscle3.8.31_i86darwin64", "-in", "sim_tempfile.fa",
        "-out", "sim_tempfile.afa", "-quiet"])

    sim_aa_seqdict = {}
    files.build_seqdict("sim_tempfile.afa",sim_aa_seqdict)
    for k in sim_aa_seqdict.keys():
        if re.search(rna_string,k):
            sim_aa_rna_seq = sim_aa_seqdict.get(k).upper()
        elif re.search(gen_string,k):
            sim_aa_gen_seq = sim_aa_seqdict.get(k).upper()
        else:
            sim_aa_ref_seq = sim_aa_seqdict.get(k).upper()

    sim_edit_list = []
    sequence.compare_aa_seqs(0,sim_aa_gen_seq,sim_aa_rna_seq,sim_aa_ref_seq,sim_edit_list)
    total_sim_score += rmath.calculate_score_diff(sim_edit_list)
    sim_scr_diffs = []
    for P,RA,GA,MA,IB,IA,SB,SA,D in sim_edit_list:
        sim_scr_diffs.append(int(D))

    if rmath.equal_vars(scr_diffs,sim_scr_diffs):
        p_val = st.ttest_ind(scr_diffs,sim_scr_diffs)
    else:
        p_val = st.ttest_ind(scr_diffs,sim_scr_diffs,equal_var=False)
    p_values.append(p_val[1])

avg_sim_score = (float(total_sim_score)/num_gens)
sig_pvals = rmath.sig_pvalue_frequency(p_values)
m_o.write("%s,%s,%.2f,%.2f,%.2f" % (name,num_edits,avg_score,avg_sim_score,sig_pvals))
m_o.write("\n")
#simulation.get_codon_percent()
#simulation.get_base_conversion()

os.remove("tempfile.fa")
os.remove("tempfile.afa")
os.remove("sim_tempfile.fa")
os.remove("sim_tempfile.afa")

# Finally, close the output file again
m_o.close()
