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
    m_o.write("name,num edits,num AA edits,average edit score,average sim score,frequency of significant edits")
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
san_ref_seq = strings.sanitize(ref_seq)
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

#print len(new_rna_seq)
#print new_gen_seq
#print new_rna_seq
num_edits,codon_positions,cumul_weights = sequence.get_positional_information(
        new_gen_seq,new_rna_seq,seq_pair.index_position())
#print num_edits
#print codon_positions
#print cumul_weights

# We don't need to use sanitized sequences, because the translate function
# removes gap characters prior to translation anyway
with open("tempfile.fa",'w') as o1:
    o1.write(">gen_seq" + "\n" + sequence.translate(gen_seq,1) + "\n")
    o1.write(">mrna_seq" + "\n" + sequence.translate(rna_seq,1) + "\n")
    o1.write(">ref_seq" + "\n" + sequence.translate(ref_seq,1) + "\n")

#print "aligning starting sequences"
subprocess.call(["muscle3.8.31_i86darwin64", "-in", "tempfile.fa",
    "-out", "tempfile.afa", "-quiet"])

aa_seqdict = {}
files.build_seqdict("tempfile.afa",aa_seqdict)
for k in aa_seqdict.keys():
    if re.search(rna_string,k):
        #print "found a rna sequence"
        aa_rna_seq = aa_seqdict.get(k).upper()
    elif re.search(gen_string,k):
        #print "found a gen sequence"
        aa_gen_seq = aa_seqdict.get(k).upper()
    else:
        #print "found a ref sequence"
        aa_ref_seq = aa_seqdict.get(k).upper()

#print rna_seq
#print
#print aa_rna_seq
#print
#print gen_seq
#print
#print aa_gen_seq

# Find the beginning and end of aligned region
i = 0
j = 0
try:
    # Compare genomic and RNA sequences to find local regions of good
    # similarity, this is taken as the start and end of aligned region
    while not sequence.compare_seqs((strings.gulp(aa_rna_seq, i, 10)),
            (strings.gulp(aa_gen_seq, i, 10)), 5):
        if i > len(aa_rna_seq):
            raise IndexError
        i += 1
    while not sequence.compare_seqs((strings.gulp(aa_rna_seq[::-1], j, 10)),
            (strings.gulp(aa_gen_seq[::-1], j, 10)), 5):
        if j > len(aa_gen_seq):
            raise IndexError
        j += 1
# If we get an index error then we cannot find start and end of both sequences
except(IndexError):
    print "Could not discern aligned part of sequences"
    # Exit cleanly
    sys.exit(0)


edit_list = []
sequence.compare_aa_seqs(0,aa_gen_seq[i:(len(aa_gen_seq)-j)],
    aa_rna_seq[i:(len(aa_rna_seq)-j)],
    aa_ref_seq[i:(len(aa_ref_seq)-j)],edit_list)
avg_score = rmath.calculate_score_diff(edit_list)
scr_diffs = []
for P,RA,GA,MA,IB,IA,SB,SA,D in edit_list:
    scr_diffs.append(int(D))
if len(scr_diffs) == 0:
    print "No amino acid changes to compare to"
    # Exit cleanly
    sys.exit(0)

p_values = []
simulation = classes.Simulation()
total_sim_score = 0
#gen_list = [b for b in new_gen_seq]
gen_list = []
for rg,rm in zip(new_gen_seq,new_rna_seq):
    # Added a case to deal with frameshifts
    if rg != '-' and rm == '-':
        gen_list.append('-')
    else:
        gen_list.append(rg)
#print gen_list
#print len(gen_list)
#print "starting simulations"
x = 0
while x < num_gens:
    new_seq = "".join(sequence.weighted_mutation(gen_list,num_edits,
        codon_positions,cumul_weights,simulation))
    #print new_seq
    #print len(new_seq)
    edited_rna_seq = rna_start + new_seq + rna_end
    #print len(rna_seq)
    #print len(edited_rna_seq)

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

    # Find the beginning and end of aligned region
    i = 0
    j = 0
    try:
        # Compare genomic and RNA sequences to find local regions of good
        # similarity, this is taken as the start and end of aligned region
        while not sequence.compare_seqs((strings.gulp(sim_aa_rna_seq, i, 10)),
            (strings.gulp(sim_aa_gen_seq, i, 10)), 5):
            if i > len(sim_aa_rna_seq):
                raise IndexError
            i += 1
        while not sequence.compare_seqs((strings.gulp(sim_aa_rna_seq[::-1], j, 10)),
            (strings.gulp(sim_aa_gen_seq[::-1], j, 10)), 5):
            if j > len(sim_aa_rna_seq):
                raise IndexError
            j += 1
    # If we get an index error then we cannot find start and end of both sequences
    except(IndexError):
        pass

    sim_edit_list = []
    sequence.compare_aa_seqs(0,sim_aa_gen_seq[i:(len(sim_aa_gen_seq)-j)],
            sim_aa_rna_seq[i:(len(sim_aa_rna_seq)-j)],
            sim_aa_ref_seq[i:(len(sim_aa_ref_seq)-j)],sim_edit_list)
    total_sim_score += rmath.calculate_score_diff(sim_edit_list)
    sim_scr_diffs = []
    for P,RA,GA,MA,IB,IA,SB,SA,D in sim_edit_list:
        sim_scr_diffs.append(int(D))
    if len(sim_scr_diffs) != 0:
        #print "calculating stats"
        #print name
        #print len(scr_diffs)
        #print len(sim_scr_diffs)
        #print
        p_val = st.ranksums(scr_diffs,sim_scr_diffs)
        p_values.append(p_val[1])
        x += 1
    else:
        #print "skipping"
        pass

num_exp_edits = len(scr_diffs)
avg_sim_score = (float(total_sim_score)/num_gens)
sig_pvals = rmath.sig_pvalue_frequency(p_values)
m_o.write("%s,%s,%s,%.2f,%.2f,%.2f" % (name,num_edits,num_exp_edits,avg_score,avg_sim_score,sig_pvals))
m_o.write("\n")
#simulation.get_codon_percent()
#simulation.get_base_conversion()

os.remove("tempfile.fa")
os.remove("tempfile.afa")
os.remove("sim_tempfile.fa")
os.remove("sim_tempfile.afa")

# Finally, close the output file again
m_o.close()
