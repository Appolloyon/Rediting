#!/usr/bin/env python

import re
import os
import sys
import argparse

from rediting.classes import classes,matrices
from rediting.util import files, rmath, sequence, strings

parser = argparse.ArgumentParser(
    description = """Calculates editing stats between genomic/RNA sequences""",
    epilog = """This program assumes that the genomic and RNA sequences for a
    given gene are provided in the same file in an aligned FASTA format, such as that
    output by MAFFT or MUSCLE. It will go through and calculate information
    regarding the editing events, as identified by differences in the aligned
    sequences. These include nucleotide, codon, and amino acid changes, as well
    as summarizing changes in all of these as well as GC content. Depending on
    user input, one or more files will be created either in .txt or .csv format.
    In order to distinguish between genomic and RNA sequences, the user must
    specify a distinguishing string (word or list of characters) present in
    the FASTA header of each (e.g. 'RNA' or 'mRNA' for RNA sequences.""")
parser.add_argument('-in', '--infile', help='infile with aligned sequences')
parser.add_argument('-out', '--outfile', help='name for master outfile')
parser.add_argument('-n', '--name', help='name to append to output file')
parser.add_argument('-g', '--gene', help='gene name for output files')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence header')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequence headers')
parser.add_argument('-neq', '--numequal', help='number of equal residues out of "size"\
        to signify start/end of alignment', default=7)
parser.add_argument('-s', '--size', help='number of residues to compare to determine start/end\
        of an alignment', default=9)
parser.add_argument('-e', '--edits', action='store_true', help='summarize editing types')
parser.add_argument('-c', '--codon', action='store_true', help='summarize codon usage difference')
parser.add_argument('-t', '--polyt', action='store_true', help='calculate polyT')
parser.add_argument('-p', '--percent', help='percent cut-off for polyT', default=70)
args = parser.parse_args()

# Store argparse variables in global variables for later use
num_equal = int(args.numequal)
size = int(args.size)
percent = int(args.percent)
name = args.name
gene = args.gene

# Create a "master" outfile to collate data from multiple files
m_out = args.outfile
m_short = m_out.rsplit('.',1)[0] # Used for other 'master' files
if os.path.isfile(m_out):
    # Appends if specified file already exists
    m_o = files.get_variable_file_handle(m_out,'a')
else:
    # First time file is opened, write header lines
    mlist = ['gene','GC before','GC after','nucleotide length','amino acid length',
        'number edits','percent edits','first position edits','second position edits',
        'third position edits','percent edits in first two positions',
        'percent non-synonymous edits','number amino acid edits','percent amino acid edits',
        'average edit score']
    if args.polyt:
        mlist.extend(['fraction polyT before','fraction polyT after',('fraction ' + str(percent) +
            ' percent polyT before'),('fraction ' + str(percent) + ' percent polyT after')])
    m_o = files.get_variable_file_handle(m_out,'w',',',mlist)

# In addition to the master file, create separate outfile(s)
# depending on which program conditions are specified
b_out = name + "_basic_editing.csv"

# Create a GC outfile to collate codon-position-specific GC data from multiple files
t_out = m_short + "_GC_position.csv"
if os.path.isfile(t_out):
    t_o = files.get_variable_file_handle(t_out,'a')
else:
    tlist = ['gene','1st GC before','1st GC after','2nd GC before',
        '2nd GC after','3rd GC before','3rd GC after']
    t_o = files.get_variable_file_handle(t_out,'w',',',tlist)

if args.edits:
    # Create files to collate positional edit data from multiple files
    e_out = m_short + "_editing_types.csv"
    p1_out = m_short + "_first_pos_editing_types.csv"
    p2_out = m_short + "_second_pos_editing_types.csv"
    p3_out = m_short + "_third_pos_editing_types.csv"
    # Same headers for each file
    plist = ['gene','A to T','A to G','A to C','T to A','T to G','T to C',
        'G to A','G to T','G to C','C to A','C to T','C to G']
    # Overall
    if os.path.isfile(e_out):
        e_o = files.get_variable_file_handle(e_out,'a')
    else:
        e_o = files.get_variable_file_handle(e_out,'w',',',plist)
    # First codon position
    if os.path.isfile(p1_out):
        p1_o = files.get_variable_file_handle(p1_out,'a')
    else:
        p1_o = files.get_variable_file_handle(p1_out,'w',',',plist)
    # Second codon position
    if os.path.isfile(p2_out):
        p2_o = files.get_variable_file_handle(p2_out,'a')
    else:
        p2_o = files.get_variable_file_handle(p2_out,'w',',',plist)
    # Third codon position
    if os.path.isfile(p3_out):
        p3_o = files.get_variable_file_handle(p3_out,'a')
    else:
        p3_o = files.get_variable_file_handle(p3_out,'w',',',plist)

if args.codon:
    c_out = name + "_codon_preference.csv"

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

# We directly compare aligned sequences, but class implementation uses
# unaligned sequences (i.e. no gap characters '-')
san_rna_seq = strings.sanitize(rna_seq)
san_gen_seq = strings.sanitize(gen_seq)
seq_pair = classes.SeqPair(san_rna_seq,san_gen_seq,name)

# Find beginning and end of aligned region
i = 0
j = 0
try:
    # Compare genomic and RNA sequences to find local regions of good
    # similarity, this is taken as the start and end of aligned region
    while not sequence.compare_seqs((strings.gulp(rna_seq, i, size)),
            (strings.gulp(gen_seq, i, size)), num_equal):
        # If we find residues in either sequence, we need to increment
        # certain class values accordingly
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
    print "Could not discern aligned part of sequences for gene " + str(name)
    # Exit cleanly
    sys.exit(0)

# Once we know the start and end, simply chop off everything else
new_rna_seq = rna_seq[i:(len(rna_seq)-j)]
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]
# The aligned amino acid length is taken as the length of the conceptual
# translation of the genomic sequence. Before it was just taken as this
# length divided by 3, but this is more accurate
aa_length = len(sequence.translate(new_gen_seq,seq_pair.index_position()))

edit_list = []
if args.edits:
    num_edited_res = 0

# Make sequences for each codon position
codon_pos_seq_dict = {
    'gen_first':'',
    'gen_second':'',
    'gen_third':'',
    'rna_first':'',
    'rna_second':'',
    'rna_third':'',
}

# Compare matching regions and look for unequal residues (edits)
for i, (rg, rm) in enumerate(zip(new_gen_seq, new_rna_seq)):
    # Regardless, we need to update the sequences
    sequence.update_pos_seq_dict(codon_pos_seq_dict,rg,rm,seq_pair.codon_pos)
    # There is an insertion in the RNA only
    if rg == '-' and rm != '-':
        seq_pair.incr_mrna()
    # There is an insertion in the DNA only
    elif rm == '-' and rg != '-':
        seq_pair.incr_all()
    # Base present in both, but no editing
    elif rg == rm:
        if args.codon:
            # Saves us from updating same codon more than once
            if seq_pair.codon_pos != 3 or i < 2:
                pass
            else:
                seq_pair.update_gcodons()
                seq_pair.update_mcodons()
        seq_pair.incr_all()
        seq_pair.incr_mrna()
    # Base present in both, but there is an edit!
    elif rg != rm:
        pos = seq_pair.index_nuc() + 1 # Index is different than position
        cpos = seq_pair.index_position()
        gnuc = seq_pair.lookup_gnuc()
        mnuc = seq_pair.lookup_mnuc()
        gcod = seq_pair.lookup_gcodon()
        mcod = seq_pair.lookup_mcodon()
        gaa = seq_pair.lookup_gaa()
        maa = seq_pair.lookup_maa()
        scr = (matrices.Blosum62(gaa, maa).sub_score())

        # We can identify whether the residue is present in a region of local
        # 'T' concentration, i.e. polyT
        if args.polyt:
            # Test whether the base is in a region of 4 or more sequential 'T's
            is_polyt = "N"
            # Only look at first seven bases at the start
            if i <= 4:
                polyt_test_seq = strings.gulp(new_gen_seq, 0, 7)
            # Only look at last seven bases at the end
            elif i >= len(new_gen_seq) - 4:
                polyt_test_seq = strings.gulp(new_gen_seq,
                        len(new_gen_seq)-7, 7)
            # In the middle take 3 bases on either side (seven total)
            else:
                polyt_test_seq = strings.gulp(new_gen_seq, i-3, 7)
            # Determine whether the region fits the definition of "polyT"
            if sequence.polyT(polyt_test_seq):
                is_polyt = "Y"

            # Test whether the base is present in region of X % 'T'
            is_polyt_percent = "N"
            percent_polyt_seqs = []
            for y in range(10): # i.e. for 10 base window
                try:
                    percent_polyt_seqs.append(
                            strings.gulp(new_gen_seq, i-y, 10))
                    # This should return a list of multiple overlapping
                    # 10 base windows to test for polyT
                # Should only fail towards end of sequence
                except(IndexError):
                    pass
            # Check whether any of the identified windows have X % 'T'
            if rmath.ispolyTpercent(percent_polyt_seqs, percent):
                is_polyt_percent = "Y"

        if not args.polyt:
            edit_list.append([pos,cpos,gnuc,mnuc,gcod,mcod,gaa,maa,scr])
        elif args.polyt:
            edit_list.append([pos,cpos,gnuc,mnuc,gcod,mcod,gaa,maa,scr,
                is_polyt,is_polyt_percent])

        # Only update codons if we care about them
        if args.codon:
            if seq_pair.codon_pos != 3 or i < 2:
                pass
            else:
                seq_pair.update_gcodons()
                seq_pair.update_mcodons()

        # Only update transitions if we care about them
        if args.edits:
            seq_pair.update_transdict()
            if seq_pair.codon_pos == 1:
                seq_pair.update_first_pos_transdict()
            elif seq_pair.codon_pos == 2:
                seq_pair.update_second_pos_transdict()
            elif seq_pair.codon_pos == 3:
                seq_pair.update_third_pos_transdict()
            num_edited_res += 1

        seq_pair.incr_all()
        seq_pair.incr_mrna()

# This next part simply counts the number of polyT regions in genomic and RNA
# sequences. Since this does not rely on presence/absence of edits per se, simplest
# thing to do is to run through each sequence again and count total occurrence
if args.polyt:
    num_gen_polyt = 0.0
    num_rna_polyt = 0.0
    # Generate all overlapping 7bp windows
    polyt_indices = sequence.get_indices(new_gen_seq, 7)

    for start,end in polyt_indices:
        if sequence.polyT(new_gen_seq[start:end]):
            num_gen_polyt += 1.0
        if sequence.polyT(new_rna_seq[start:end]):
            num_rna_polyt += 1.0

    num_gen_percent_polyt = 0.0
    num_rna_percent_polyt = 0.0
    # Generate all overlapping 10bp windows
    percent_polyt_indices = sequence.get_indices(new_gen_seq, 10)

    for start,end in percent_polyt_indices:
        if sequence.polyTpercent(new_gen_seq[start:end], percent):
            num_gen_percent_polyt += 1.0
        if sequence.polyTpercent(new_rna_seq[start:end], percent):
            num_rna_percent_polyt += 1.0

        # Calculate the fraction, since we take overlapping indices
        # Specifically this is the fraction of all possible windows
        # and so is related to, but distinct from, sequence length
        fraction_gen_polyt = (num_gen_polyt/len(polyt_indices)) * 100
        fraction_rna_polyt = (num_rna_polyt/len(polyt_indices)) * 100
        fraction_gen_percent_polyt = (num_gen_percent_polyt/
                len(percent_polyt_indices)) * 100
        fraction_rna_percent_polyt = (num_rna_percent_polyt/
                len(percent_polyt_indices)) * 100

subscore = 0
num_aaedits = 0
num_first_pos = 0
num_second_pos = 0
num_third_pos = 0
with open(b_out,'w') as b_o:
    # Write header line for each gene sheet
    b_o.write("position,codon position,genome base,mRNAbase,genome codon,"
        "mRNA codon,genome amino acid,mRNA amino acid,substitution score")
    if args.polyt:
        # Add columns if polyT
        b_o.write(",in a polyT tract,in a tract with " +
                str(percent) + " percent T residues")
    b_o.write("\n" * 2)

    if args.polyt:
        for P, C, GN, MN, GC, MC, GA, MA, S, IP, IPP in edit_list:
            b_o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,%s"
                    % (P,C,GN,MN,GC,MC,GA,MA,S,IP,IPP) + "\n")
            # Keep a running tally of total sub score
            subscore += int(S)
            # Count certain global values
            if GA != MA:
                num_aaedits += 1
            if C == 1:
                num_first_pos += 1
            if C == 2:
                num_second_pos += 1
            if C == 3:
                num_third_pos += 1
    else:
        # Fewer values in each list entry without polyT
        for P, C, GN, MN, GC, MC, GA, MA, S in edit_list:
            b_o.write("%s,%s,%s,%s,%s,%s,%s,%s,%s" %
                    (P,C,GN,MN,GC,MC,GA,MA,S) + "\n")
            subscore += int(S)
            if GA != MA:
                num_aaedits += 1
            if C == 1:
                num_first_pos += 1
            if C == 2:
                num_second_pos += 1
            if C == 3:
                num_third_pos += 1

# Simple % GC calculation for genomic and RNA sequences
gcb = sequence.calc_gc(new_gen_seq)
gca = sequence.calc_gc(new_rna_seq)
# This is technically the length of the aligned sequences
seqlength = len(new_rna_seq)
# Amino acid length implementation changed, see line 120
#aalength = seqlength/3
numedits = float(len(edit_list))
# Percent edits as compared to the aligned region
seqedits = (numedits/seqlength) * 100
aaedits = (float(num_aaedits)/aa_length) * 100

# If no residues are edited then this will throw an error
try:
    non_syn = (float(num_aaedits)/numedits) * 100
except(ZeroDivisionError):
    non_syn = 0
# Again, same error will occur if numedits is zero

try:
    editscore = subscore/numedits
except(ZeroDivisionError):
    editscore = 0

# Again, same error will occur if numedits is zero
try:
    first_two_pos = ((num_first_pos + num_second_pos)/numedits) * 100
except(ZeroDivisionError):
    first_two_pos = 0

m_o.write("%s,%.2f,%.2f,%s,%s,%s,%.2f,%s,%s,%s,%.2f,%.2f,%s,%.2f,%.2f" % (gene,\
        gcb,gca,seqlength,aa_length,numedits,seqedits,num_first_pos,num_second_pos,\
        num_third_pos,first_two_pos,non_syn,num_aaedits,aaedits,editscore))
if args.polyt:
    m_o.write(",%.2f,%.2f,%.2f,%.2f" % (fraction_gen_polyt,fraction_rna_polyt,\
        fraction_gen_percent_polyt,fraction_rna_percent_polyt))
m_o.write("\n")

t_o.write("%s,%.2f,%.2f,%.2f,%.2f,%.2f,%.2f" % (gene,\
        sequence.calc_gc(codon_pos_seq_dict.get('gen_first')),\
        sequence.calc_gc(codon_pos_seq_dict.get('rna_first')),\
        sequence.calc_gc(codon_pos_seq_dict.get('gen_second')),\
        sequence.calc_gc(codon_pos_seq_dict.get('rna_second')),\
        sequence.calc_gc(codon_pos_seq_dict.get('gen_third')),\
        sequence.calc_gc(codon_pos_seq_dict.get('rna_third'))))
t_o.write("\n")

if args.codon:
    with open(c_out,'w') as c_o:
        c_o.write("amino acid,codon,genome usage,mRNA usage")
        c_o.write("\n")
        # Loop through both dictionaries - if we don't sort these
        # beforehand we won't necessarily get the same value for each
        # as dictionaries are not ordered
        for k1, k2 in zip(sorted(seq_pair.gnuc_aa_dict),\
                sorted(seq_pair.mnuc_aa_dict)):
            c_o.write(k1)
            for k3, k4 in zip(sorted(seq_pair.gnuc_aa_dict[k1].keys()),\
                    sorted(seq_pair.mnuc_aa_dict[k2].keys())):
                c_o.write(',' + k3 + ',' + str(seq_pair.gnuc_aa_dict[k1][k3])\
                        + ',' + str(seq_pair.mnuc_aa_dict[k2][k4]) + "\n")
            c_o.write("\n")

if args.edits:
    files.write_transition_dict(gene,seq_pair.transition_dict,e_o)
    files.write_transition_dict(gene,seq_pair.first_pos_transdict,p1_o)
    files.write_transition_dict(gene,seq_pair.second_pos_transdict,p2_o)
    files.write_transition_dict(gene,seq_pair.third_pos_transdict,p3_o)

m_o.close()
t_o.close()
if args.edits:
    e_o.close()
    p1_o.close()
    p2_o.close()
    p3_o.close()
