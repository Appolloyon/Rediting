#!/usr/bin/env python

import re
import os
import sys
import argparse
import matplotlib.pyplot as plt

from classes import classes, matrices
from util import sequence_alignment, files, rmath, sequence, strings

parser = argparse.ArgumentParser(
    description = """Compares editing frequency between aligned genomic/RNA
        sequences and a reference""",
    epilog = """This program assumes that the genomic and RNA sequence for a
    given gene are provided in the same file in an aligned FASTA format, such
    as that output by MAFFT or MUSCLE. In addition, a reference sequence in the
    same aligned file will be used for comparison. Based on a given sliding
    window size, the program will compare the genomic sequence to the RNA sequence
    and the reference in each window frame and calculate the % editing and %
    identity, respectively. In addition, a pearson correlation will be calculated
    to determine if the two trends are related along the entire length of the
    sequence (all observed windows are treated as individual data points). In order
    to distinguish between genomic and RNA sequences, the user must specify a
    distinguishing string (word or list of characters) present in the FASTA
    header of each (e.g. 'RNA' or 'mRNA' for RNA sequences.""")
parser.add_argument('-in', '--infile', help='infile with aligned sequences')
parser.add_argument('-out', '--outfile', help='name for master outfile')
parser.add_argument('-n', '--name', help='name to append to output file')
parser.add_argument('-g', '--gene', help='gene name for output files')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence header')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequence header')
parser.add_argument('-neq', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignment', default=7)
parser.add_argument('-s', '--size', help='number of residues to compare to determine\
        start/end of an alignment', default=9)
parser.add_argument('-w', '--window_size', help='size of sliding window', default=60)
parser.add_argument('-p', '--protein', action='store_true', help='compare conceptual translations as well')
parser.add_argument('-o', '--synonymous', action='store_true', help='only record synonymous edits')
args = parser.parse_args()

window_size = float(args.window_size)
num_equal = int(args.numequal)
size = int(args.size)
name = args.name
gene = args.gene

#print args.infile

m_out = args.outfile
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    m_o = open(m_out,'w')
    m_o.write("name,percent above average edits,number obs,DF (N-2),"
        "nucleotide pearson correlation,nucleotide t value")
    if args.protein:
        m_o.write(",amino acid pearson correlation,amino acid t value")
    m_o.write("\n" * 2)

seqdict = {}
files.build_seqdict(args.infile,seqdict)

rna_string = str(args.RNA)
gen_string = str(args.genomic)
for k in seqdict.keys():
    if re.search(rna_string,k):
        rna_seq = seqdict.get(k)
    elif re.search(gen_string,k):
        gen_seq = seqdict.get(k)
    else:
        ref_seq = seqdict.get(k)

if args.protein:
    san_gen_seq = strings.sanitize(gen_seq)
    san_ref_seq = strings.sanitize(ref_seq)
    ref_pair = classes.RefPair(san_ref_seq,san_gen_seq,name)
if args.synonymous:
    san_gen_seq = strings.sanitize(gen_seq)
    san_rna_seq = strings.sanitize(rna_seq)
    seq_pair = classes.SeqPair(san_rna_seq,san_gen_seq,name)

# Need to find beginning and end of aligned region
i = 0
j = 0
try:
    while not sequence.compare_seqs((strings.gulp(rna_seq, i, size)),
            (strings.gulp(gen_seq, i, size)), num_equal):
        if gen_seq[i] != '-':
            if args.protein:
                ref_pair.incr_all_gen()
            if args.synonymous:
                seq_pair.incr_all()
        if ref_seq[i] != '-':
            if args.protein:
                ref_pair.incr_all_ref()
        if rna_seq[i] != '-':
            if args.synonymous:
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

new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]
new_ref_seq = ref_seq[i:(len(ref_seq)-j)]

# Cannot do sliding window calculations if sequences are smaller than chosen
# window size; causes numerous errors in later parts of program
if (len(new_rna_seq) < window_size or
    len(new_gen_seq) < window_size or
    len(new_ref_seq) < window_size):
    print "One or more sequences are shorter than chosen window size."
    print "Please choose a smaller window size or check sequences."
    # Exit cleanly
    sys.exit(0)

# calculates % edits between genomic and RNA sequences
compstr1 = ''
for i, (rg, rm) in enumerate(zip(new_gen_seq, new_rna_seq)):
    # If there is an insertion in the RNA only
    if rg == '-' and rm != '-':
        compstr1 += str(0)
        if args.synonymous:
            seq_pair.incr_mrna()
    # If there is an insertion in the DNA only
    elif rm == '-' and rg != '-':
        compstr1 += str(0)
        if args.synonymous:
            seq_pair.incr_all()
    # Unlike editing pipeline program, there is a possibility to
    # have a gap in both sequences
    elif rg == '-' and rm == '-':
        compstr1 += str(0)
    # Base present in both but no editing
    elif rg == rm:
        compstr1 += str(0)
        if args.synonymous:
            seq_pair.incr_all()
            seq_pair.incr_mrna()
    # Base present in both, but there is an edit!
    elif rg != rm:
        if args.synonymous:
            if seq_pair.lookup_gaa() == seq_pair.lookup_maa():
                # We ignore synonymous edits
                compstr1 += str(0)
            elif seq_pair.lookup_gaa() != seq_pair.lookup_maa():
                compstr1 += str(1)
            seq_pair.incr_all()
            seq_pair.incr_mrna()
        else:
            compstr1 += str(1)
    else:
        pass

# Calculates percent edits for each window
edit_list = []
for start,end in sequence.get_indices(compstr1, window_size):
    try:
        edit_list.append(sequence.calc_percent(
            compstr1, start, end, window_size))
    except(ValueError,IndexError):
        pass

# Calculates % sequence identity beween genomic and reference sequences
compstr2 = ''
for i, (res1, res2) in enumerate(zip(new_gen_seq, new_ref_seq)):
    if res1 == res2:
        compstr2 += str(1)
    elif res1 != res2:
        compstr2 += str(0)
    else:
        pass

# Calculates percent gen/ref identity for each window
identity_list = []
for start,end in sequence.get_indices(compstr2, window_size):
    try:
        identity_list.append(sequence.calc_percent(
            compstr2, start, end, window_size))
    except(ValueError,IndexError):
        pass

if args.protein:
    similarity_list = []
    for i, (rg,rr) in enumerate(zip(new_gen_seq, new_ref_seq)):
        similarity_sum = 0.0
        if new_gen_seq[i] != '-':
            if i > 0:
                ref_pair.incr_all_gen()
        if new_ref_seq[i] != '-':
            if i > 0:
                ref_pair.incr_all_ref()

        rpos = ref_pair.index_rposition()
        gpos = ref_pair.index_gposition()
        rnuc_seq = strings.gulp(new_ref_seq, i, int(window_size))
        raa_seq = sequence.translate(rnuc_seq,rpos)
        gnuc_seq = strings.gulp(new_gen_seq, i, int(window_size))
        gaa_seq = sequence.translate(gnuc_seq,gpos)

        if (len(rnuc_seq) == int(window_size) and
                len(gnuc_seq) == int(window_size)):
            # If the lengths aren't equal, align them
            if len(raa_seq) != len(gaa_seq):
                raa_seq,gaa_seq = sequence_alignment.affine_align(raa_seq,gaa_seq)
            # Determine how similar raa_seq and gaa_seq are
            for raa,gaa in zip(raa_seq,gaa_seq):
                # Gaps are neutral
                if raa == '-' or gaa == '-':
                    similarity_sum += 0.0
                # Identities score 1
                elif raa == gaa:
                    similarity_sum += 1.0
                # Similar residues have positive scores
                elif matrices.Blosum62(raa,gaa).sub_score() > 0:
                    similarity_sum += 0.5
                else:
                    similarity_sum += 0.0

            if similarity_sum == 0.0 or len(raa_seq) == 0:
                similarity_list.append(0.0)
            else:
                try:
                    similarity_list.append(
                            (similarity_sum/float(len(raa_seq))*100))
                except(ValueError,IndexError,ZeroDivisionError):
                    pass

# Get the average of all edits over windows
try:
    edit_mean = rmath.calc_mean(edit_list)
    average_list = []
    for i in range(len(edit_list)):
        average_list.append(edit_mean)
# If no edits occurred, then edit_list is empty
# and calc_mean will throw a ZeroDivError
except(ZeroDivisionError):
    for i in range(len(edit_list)):
        average_list.append(0.0)

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
    percent_above_average_edits = 0.0
num_obs = len(edit_list)

identity_mean = rmath.calc_mean(identity_list)
try:
    PC = rmath.calc_pearson(edit_list,
            identity_list, edit_mean, identity_mean)
# This calculation can fail, especially if edit_list is 0 length
# In that case, just call it 0
except:
    PC = 0.0
tvalue = rmath.calc_tvalue(PC, num_obs)

if args.protein:
    similarity_mean = rmath.calc_mean(similarity_list)
    try:
        aa_PC = rmath.calc_pearson(edit_list,
            similarity_list, edit_mean, similarity_mean)
    except:
        aa_PC = 0.0
    aa_tvalue = rmath.calc_tvalue(aa_PC, num_obs)

#if args.long_output:
m_o.write("%s,%.2f,%d,%d,%f,%f" % (name,percent_above_average_edits,
    num_obs,(num_obs - 2),PC,tvalue))
if args.protein:
    m_o.write(",%f,%f" % (aa_PC,aa_tvalue))
m_o.write('\n')

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

ax1.plot([e for e in edit_list], color='b')
ax1.plot([mean for mean in average_list], linestyle='--', color='k')
ax2.plot([i for i in identity_list], color='m')
if args.protein:
    ax2.plot([i for i in similarity_list], color='g')

plt.title('%s' % (name))
ax1.set_xlabel('sliding window position')
ax1.set_ylabel('percent edited residues in RNA', color='b')
ax2.set_ylabel('percent sequence identity to reference', color='m')
if args.protein:
    # Extend plot to the right slightly
    fig.subplots_adjust(right=0.88)
    # Add some additional text to the ylabel
    plt.figtext(0.95,0.77,'percent amino acid similarity to reference', color='g',
        rotation='vertical')

if args.protein:
    textstr = ("percent above average = %.2f\nnumber of windows = %d\n"
        "nucleotide pearson's correlation = %.2f\nnucleotide tvalue = %.2f\n"
        "amino acid pearson's correlation = %.2f\namino acid tvalue = %.2f" %
        (percent_above_average_edits,num_obs,PC,tvalue,aa_PC,aa_tvalue))
else:
    textstr = ("percent above average = %.2f\nnumber of windows = %d\n"
        "nucleotide pearson's correlation = %.2f\nnucleotide tvalue = %2.f" %
        (percent_above_average_edits,PC,num_obs,tvalue))
props = dict(boxstyle='round', facecolor='white', alpha=0.75)
ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes,
        fontsize=8, verticalalignment='top', bbox=props)

plt.savefig('%s.pdf' % (name))
plt.close()

#if args.long_output:
m_o.close()
