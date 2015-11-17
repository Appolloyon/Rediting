#!/usr/bin/env python

import re
import os
import argparse
import matplotlib.pyplot as plt

from classes import RefPair
from matrices import Blosum62
from sequence_alignment import affine_align
from functions import gulp, compare_seqs, sanitize, build_seqdict, calc_percent,\
        get_indices, calc_mean, calc_pearson, calc_tvalue, translate

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
args = parser.parse_args()

window_size = float(args.window_size)
num_equal = int(args.numequal)
size = int(args.size)
name = args.name
gene = args.gene

print args.infile

m_out = args.outfile
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    m_o = open(m_out,'w')
    m_o.write("name,percent above average edits,number obs,DF (N-2),nucleotide pearson correlation,nucleotide t value")
    if args.protein:
        m_o.write(",amino acid pearson correlation,amino acid t value")
    m_o.write("\n" * 2)

#for infile in args.infiles:
    #name = infile.split('.')[0]
    #gene = name.split('_')[1]

seqdict = {}
build_seqdict(args.infile,seqdict)

rna_string = str(args.RNA)
gen_string = str(args.genomic)
for k in seqdict.keys():
    if re.search(rna_string,k):
        rna_seq = seqdict.get(k)
    elif re.search(gen_string,k):
        gen_seq = seqdict.get(k)
    else:
        ref_seq = seqdict.get(k)
#print rna_seq
#print gen_seq

if args.protein:
    san_gen_seq = sanitize(gen_seq)
    san_ref_seq = sanitize(ref_seq)
    ref_pair = RefPair(san_ref_seq,san_gen_seq,name)

"""
Another insidious bug, this comparison can effectively go on
past the length of the sequence and throw an IndexError
"""
i = 0
j = 0
while not compare_seqs((gulp(rna_seq, i, size)),
        (gulp(gen_seq, i, size)), num_equal):  #start of alignment
    #print i
    #print gen_seq[i]
    #print ref_seq[i]
    if gen_seq[i] != '-':
        ref_pair.incr_all_gen()
    if ref_seq[i] != '-':
        ref_pair.incr_all_ref()
    i += 1
while not compare_seqs((gulp(rna_seq[::-1], j, size)),
        (gulp(gen_seq[::-1], j, size)), num_equal):  #end of alignment
    #print j
    j += 1

new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
print rna_seq
print new_rna_seq
print
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]
print gen_seq
print new_gen_seq
print
new_ref_seq = ref_seq[i:(len(ref_seq)-j)]
print ref_seq
print new_ref_seq
print len(new_ref_seq)
print

# calculates % edits between genomic and RNA sequences
compstr1 = ''
for i, (res1, res2) in enumerate(zip(new_rna_seq, new_gen_seq)):
    if (res1 == '-' or res2 == '-') or res1 == res2:
        compstr1 += str(0)
    elif res1 != res2:
        compstr1 += str(1)
    else:
        pass
#print compstr1

"""
Noted a very insidious bug in the program, namely that if the aligned sequence
length is less than the window size, there is nothing to compare and the
program will throw a ZeroDivisionError later on when it attemps to calulcate the
mean of an empty list
"""
edit_list = []
for start,end in get_indices(compstr1, window_size):
    try:
        edit_list.append(calc_percent(compstr1, start, end, window_size))
    except(ValueError,IndexError):
        pass
#print edit_list

# calculates % sequence identity beween genomic and reference sequences
compstr2 = ''
for i, (res1, res2) in enumerate(zip(new_gen_seq, new_ref_seq)):
    if res1 == res2:
        compstr2 += str(1)
    elif res1 != res2:
        compstr2 += str(0)
    else:
        pass
#print compstr2

identity_list = []
for start,end in get_indices(compstr2, window_size):
    try:
        identity_list.append(calc_percent(compstr2, start, end, window_size))
    except(ValueError,IndexError):
        pass
#print identity_list

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
        rnuc_seq = gulp(new_ref_seq, i, int(window_size))
        #print "reference nucleotide sequence:" + rnuc_seq
        raa_seq = translate(rnuc_seq,rpos)
        #print "reference amino acid sequence:" + raa_seq
        gnuc_seq = gulp(new_gen_seq, i, int(window_size))
        #print "genomic nucleotide sequence:" + gnuc_seq
        gaa_seq = translate(gnuc_seq,gpos)
        #print "genomic amino acid sequence:" + gaa_seq

        if len(rnuc_seq) == int(window_size) and len(gnuc_seq) == int(window_size):
            if len(raa_seq) != len(gaa_seq): # sequence require further alignment
                raa_seq,gaa_seq = affine_align(raa_seq,gaa_seq)
                #print raa_seq
                #print gaa_seq
                #print
            for raa,gaa in zip(raa_seq,gaa_seq):
                if raa == '-' or gaa == '-':
                    similarity_sum += 0.0
                elif raa == gaa:
                    similarity_sum += 1.0
                elif Blosum62(raa,gaa).sub_score() > 0:
                    similarity_sum += 0.5
                else:
                    similarity_sum += 0.0

            if similarity_sum == 0.0 or len(raa_seq) == 0:
                similarity_list.append(0.0)
            else:
                try:
                    similarity_list.append((similarity_sum/float(len(raa_seq))*100))
                except(ValueError,IndexError,ZeroDivisionError):
                    pass

try:
    edit_mean = calc_mean(edit_list)
    average_list = []
    for i in range(len(edit_list)):
        average_list.append(edit_mean)
except(ZeroDivisionError):
    for i in range(len(edit_list)):
        average_list.append(0.0)

edits_above_average = 0.0
total_edits = 0.0
for edit in edit_list:
    if edit > edit_mean:
        total_edits += edit
        edits_above_average += edit
    else:
        total_edits += edit

try:
    percent_above_average_edits = (edits_above_average/total_edits) * 100
except:
    percent_above_average_edits = 0.0
num_obs = len(edit_list)

print identity_list
identity_mean = calc_mean(identity_list)
try:
    PC = calc_pearson(edit_list, identity_list, edit_mean, identity_mean)
except:
    PC = 0.0
tvalue = calc_tvalue(PC, num_obs)

if args.protein:
    similarity_mean = calc_mean(similarity_list)
    aa_PC = calc_pearson(edit_list, similarity_list, edit_mean, similarity_mean)
    aa_tvalue = calc_tvalue(aa_PC, num_obs)

#if args.long_output:
m_o.write("%s,%.2f,%d,%d,%f,%f" % (name,percent_above_average_edits,num_obs,(num_obs - 2),PC,tvalue))
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
    textstr = "percent above average = %.2f\nnumber of windows = %d\nnucleotide pearson's correlation = %.2f\nnucleotide tvalue = %.2f\namino acid pearson's correlation = %.2f\namino acid tvalue = %.2f" % (percent_above_average_edits,num_obs,PC,tvalue,aa_PC,aa_tvalue)
else:
    textstr = "percent above average = %.2f\nnumber of windows = %d\nnucleotide pearson's correlation = %.2f\nnucleotide tvalue = %2.f" % (percent_above_average_edits,PC,num_obs,tvalue)
props = dict(boxstyle='round', facecolor='white', alpha=0.75)
ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=8, verticalalignment='top', bbox=props)

plt.savefig('%s.pdf' % (name))
plt.close()

#if args.long_output:
m_o.close()
