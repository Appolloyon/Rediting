#!/usr/bin/env python

"""
This program reads in aligned files with mRNA and DNA sequences of a gene from
one organism and the corresponding genome sequence from Emiliania. It calculates
the percent edits between the transcript and the DNA sequence, and then
compares the sequence identity between the two genome sequences over the same
range. Both of these values are graphed and a pearson's correlation is calculated.

Changelog:
Created June 25, 2015
"""

import re
import argparse
import matplotlib.pyplot as plt

from functions import gulp, compare_seqs, build_seqdict, calc_percent,\
        get_indices, calc_mean, calc_pearson, calc_tvalue

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
parser.add_argument('infiles', nargs='+', help='list of infiles')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence header')
parser.add_argument('-g', '--genomic', help='unique string present in genomic sequence header')
parser.add_argument('-n', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignment', default=7)
parser.add_argument('-s', '--size', help='number of residues to compare to determine\
        start/end of an alignment', default=9)
parser.add_argument('-w', '--window_size', help='size of sliding window', default=60)
parser.add_argument('-o', '--long_output', action='store_true', help='output summary csv file')
args = parser.parse_args()

window_size = float(args.window_size)

if args.long_output:
    m_out = "master_sliding_window_out.csv"
    m_o = open(m_out, 'w')
    m_o.write("name,percent above average edits,number obs,DF (N-2),pearson correlation,t value")
    m_o.write("\n" * 2)

for infile in args.infiles:
    #name = ((os.path.basename((file)).strip(".afa")))
    name = infile.split('.')[0]
    gene = name.split('_')[1]

    seqdict = {}
    build_seqdict(infile,seqdict)

    rna_string = str(args.RNA)
    gen_string = str(args.genomic)
    for k in seqdict.keys():
        if re.search(rna_string,k):
            rna_seq = seqdict.get(k)
        elif re.search(gen_string,k):
            gen_seq = seqdict.get(k)
        else:
            ref_seq = seqdict.get(k)

    num_equal = int(args.numequal)
    size = int(args.size)
    i = 0
    j = 0
    while not compare_seqs((gulp(rna_seq, i, size)),
            (gulp(gen_seq, i, size)), num_equal):  #start of alignment
        i += 1
    while not compare_seqs((gulp(rna_seq[::-1], j, size)),
            (gulp(gen_seq[::-1], j, size)), num_equal):  #end of alignment
        j += 1

    new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
    new_gen_seq = gen_seq[i:(len(gen_seq)-j)]
    new_ref_seq = ref_seq[i:(len(ref_seq)-j)]

    # calculates % edits between genomic and RNA sequences
    compstr1 = ''
    for i, (res1, res2) in enumerate(zip(new_rna_seq, new_gen_seq)):
        if (res1 == '-' or res2 == '-') or res1 == res2:
            compstr1 += str(0)
        elif res1 != res2:
            compstr1 += str(1)
        else:
            pass

    edit_list = []
    for start,end in get_indices(compstr1, window_size):
        try:
            edit_list.append(calc_percent(compstr1, start, end, window_size))
        except(ValueError,IndexError):
            pass

    # calculates % sequence identity beween genomic and reference sequences
    compstr2 = ''
    for i, (res1, res2) in enumerate(zip(new_gen_seq, new_ref_seq)):
        if res1 == res2:
            compstr2 += str(1)
        elif res1 != res2:
            compstr2 += str(0)
        else:
            pass

    identity_list = []
    for start,end in get_indices(compstr2, window_size):
        try:
            identity_list.append(calc_percent(compstr2, start, end, window_size))
        except(ValueError,IndexError):
            pass

    edit_mean = calc_mean(edit_list)
    average_list = []
    for i in range(len(edit_list)):
        average_list.append(edit_mean)
    identity_mean = calc_mean(identity_list)

    edits_above_average = 0.0
    total_edits = 0.0
    for edit in edit_list:
        if edit > edit_mean:
            total_edits += edit
            edits_above_average += edit
        else:
            total_edits += edit

    percent_above_average_edits = (edits_above_average/total_edits) * 100
    #print percent_above_average_edits
    PC = calc_pearson(edit_list, identity_list, edit_mean, identity_mean)
    #print PC
    num_obs = len(edit_list)
    #print num_obs
    tvalue = calc_tvalue(PC, num_obs)
    #print tvalue

    if args.long_output:
        m_o.write("%s,%.2f,%d,%d,%f,%f" % (name,percent_above_average_edits,num_obs,(num_obs - 2),PC,tvalue) + "\n")

    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()

    ax1.plot([e for e in edit_list], color='b')
    ax1.plot([mean for mean in average_list], linestyle='--', color='k')
    ax2.plot([i for i in identity_list], color='m')

    plt.title('%s' % (name))
    ax1.set_xlabel('sliding window position')
    ax1.set_ylabel('percent edited residues in RNA', color='b')
    ax2.set_ylabel('percent sequence identity to reference', color='m')

    textstr = "percent above average = %.2f\npearson's correlation = %.2f\nnumber of windows = %d\ntvalue = %2.f"%(percent_above_average_edits,PC,num_obs,tvalue)
    props = dict(boxstyle='round', facecolor='white', alpha=0.75)
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=8, verticalalignment='top', bbox=props)

    plt.savefig('%s.pdf' % (name))
    plt.close()

if args.long_output:
    m_o.close()
