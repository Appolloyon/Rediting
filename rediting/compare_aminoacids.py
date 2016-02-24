#!/usr/bin/env python

import os
import re
import sys
import argparse

from classes import matrices
from util import files, sequence, strings

parser = argparse.ArgumentParser(
    description = """Calculates amino acids before and after editing to a reference""",
    epilog = """Write some kind of description here.""")
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

m_out = args.outfile
# appends if specified file already exists
if os.path.isfile(m_out):
    m_o = open(m_out,'a')
else:
    m_o = open(m_out,'w')
    m_o.write("gene,num AA changes,num identical before,num identical after,"
        "num similar before,num similar after,avgerage edit score diff")
    m_o.write("\n" * 2)

b_out = name + "_aminoacid_changes.csv"

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
    while not sequence.compare_seqs((strings.gulp(rna_seq, i, size)),
            (strings.gulp(gen_seq, i, size)), num_equal):
        if gen_seq[i] != '-':
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

new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
#print new_rna_seq
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]
#print new_gen_seq
new_ref_seq = ref_seq[i:(len(ref_seq)-j)]

edit_list = []

# calculates amino acid edits
for i, (rg,rm) in enumerate(zip(new_gen_seq, new_rna_seq)):
    ident_before = 'N'
    ident_after = 'N'
    sim_before = 'N'
    sim_after = 'N'
    # insert in mrna sequence
    if rg == '-' and rm != '-':
        pass
    # insert in genomic sequence
    elif rg != '-' and rm == '-':
        # don't do anything, but move the index along
        gen_index += 1
    # gaps in both sequences, no edits possible
    elif rg == '-' and rm == '-':
        pass
    # no edits, but move the index along
    elif rg == rm:
        gen_index += 1
    # edit detected
    elif rg != rm:
        ref_aa = new_ref_seq[i]
        # we only care about edits we can actually compare
        if ref_aa != '-':
            if rg == ref_aa:
                ident_before = 'Y'
            if rm == ref_aa:
                ident_after = 'Y'
            scr_before = matrices.Blosum62(rg,ref_aa).sub_score()
            scr_after = matrices.Blosum62(rm,ref_aa).sub_score()
            if scr_before > 0:
                sim_before = 'Y'
            if scr_after > 0:
                sim_after = 'Y'
            scr_diff = (scr_after - scr_before)

            edit_list.append([(gen_index+1),ref_aa,rg,rm,ident_before,
                ident_after,sim_before,sim_after,scr_diff])

        gen_index += 1

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
        diff_score += int(D)
        if IB == 'Y':
            num_ident_before += 1
        if IA == 'Y':
            num_ident_after += 1
        if SB == 'Y':
            num_similar_before += 1
        if SA == 'Y':
            num_similar_after += 1

numedits = float(len(edit_list))
try:
    avg_diff = diff_score/numedits
except(ZeroDivisionError):
    avg_diff = 0

m_o.write("%s,%s,%s,%s,%s,%s,%.2f" % (gene,numedits,num_ident_before,\
        num_ident_after,num_similar_before,num_similar_after,avg_diff))
m_o.write("\n")

m_o.close()
