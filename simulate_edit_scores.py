#!/usr/bin/env python

import os
import re
import sys
import argparse
import subprocess
import scipy.stats as st

from rediting.classes import classes
from rediting.util import files, sequence, strings, rmath

parser = argparse.ArgumentParser(
    description = """Simulates an equivalent number of edits with similar
    edit type and codon position preferences as real sequences. Then compares
    conceptual translations to determine whether editing is more corrective
    than expected by chance.""",
    epilog = """This program requires a pre-aligned file with genomic, RNA, and
    reference sequences. It then calculates all edits in the aligned region of
    interest between all three sequences, and simulates editing events based on
    the observed preferences. Finally, the sequences are translated and the
    edit scores calculated, as performed in the compare_aminoacids.py script
    from this same package. Lastly, for each simulation, these scores are
    compared using a ranksums test to determine significance.""")
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

# Create a "master" outfile to collate data from multiple runs
m_out = args.outfile
s_out = m_out.rsplit('.',1)[0] + "_stats.csv"
if os.path.isfile(m_out):
    # Appends if specified file already exists
    m_o = files.get_variable_file_handle(m_out,'a')
else:
    # The first time the file is opened, write header lines
    mlist = ['gene','number nucleotide edits','number AA edits','average number sim AA edits',
            'average edit score','average sim edit score','frequency of significant editing']
    m_o = files.get_variable_file_handle(m_out,'w',',',mlist)
if os.path.isfile(s_out):
    s_o = files.get_variable_file_handle(s_out,'a')
else:
    slist = ['gene','num 1st pos','num 2nd pos','num 3rd pos','A to T','A to G',
            'A to C','T to A','T to G','T to C','G to A','G to T','G to C',
            'C to A','C to T','C to G']
    s_o = files.get_variable_file_handle(s_out,'w',',',slist)

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

# We directly compare aligned sequences, but class implementation uses
# unaligned sequences (i.e. no gap characters '-')
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
# Save the first and last parts of the RNA sequence to glue back on
# later once the sequence has been "edited"
rna_start = rna_seq[:i]
rna_end = rna_seq[(len(rna_seq)-j):]
new_gen_seq = gen_seq[i:(len(gen_seq)-j)]

# In order to "edit" the sequence accurately, need to know the distribution
# of edits among codon positions and base conversions
num_edits,codon_positions,cumul_weights = sequence.get_positional_information(
        new_gen_seq,new_rna_seq,seq_pair.index_position())

#print num_edits
#print codon_positions
for weight in cumul_weights:
    print weight

# In order for MUSCLE to align, write sequences to a temporary file
# We don't need to use sanitized sequences, because the translate function
# removes gap characters prior to translation anyway
with open("tempfile.fa",'w') as o1:
    o1.write(">gen_seq" + "\n" + sequence.translate(gen_seq,1) + "\n")
    o1.write(">mrna_seq" + "\n" + sequence.translate(rna_seq,1) + "\n")
    o1.write(">ref_seq" + "\n" + sequence.translate(ref_seq,1) + "\n")

# External call to align the file with MUSCLE
subprocess.call(["muscle3.8.31_i86darwin64", "-in", "tempfile.fa",
    "-out", "tempfile.afa", "-quiet"])

# Now these sequences are amino acid translations
aa_seqdict = {}
files.build_seqdict("tempfile.afa",aa_seqdict)
for k in aa_seqdict.keys():
    if re.search(rna_string,k):
        aa_rna_seq = aa_seqdict.get(k).upper()
    elif re.search(gen_string,k):
        aa_gen_seq = aa_seqdict.get(k).upper()
    else:
        aa_ref_seq = aa_seqdict.get(k).upper()

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
# Determine the consequence of editing at the amino acid level
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

# Initialize important variables for the simulation
p_values = []
simulation = classes.Simulation()
total_sim_score = 0
total_sim_aa_edits = 0
gen_list = []
# In order to change bases in place, need to convert the sequence to a list
# This used to be a simple list comp expression but frameshifts necessitate
# comparing both genomic and rna sequences
for rg,rm in zip(new_gen_seq,new_rna_seq):
    # Added a case to deal with frameshifts
    if rg != '-' and rm == '-':
        gen_list.append('-')
    # Otherwise, just append
    else:
        gen_list.append(rg)
x = 0
while x < num_gens:
    # "Edit" the genomic sequence over the aligned region by applying edits to
    # it as per the patterns observed for the real sequences
    new_seq = "".join(sequence.weighted_mutation(gen_list,num_edits,
        codon_positions,cumul_weights,simulation))
    edited_rna_seq = rna_start + new_seq + rna_end

    # Same as above, need to align the sequences so we write to a file
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
    total_sim_aa_edits += len(sim_edit_list)
    sim_scr_diffs = []
    for P,RA,GA,MA,IB,IA,SB,SA,D in sim_edit_list:
        sim_scr_diffs.append(int(D))
    # If there are no edits then we need to try again
    if len(sim_scr_diffs) != 0:
        p_val = st.ranksums(scr_diffs,sim_scr_diffs)
        p_values.append(p_val[1])
        # By incrementing the counter we are saying that we
        # obtained edits for comparison
        x += 1
    else:
        #print "skipping"
        pass

num_exp_edits = len(scr_diffs)
avg_sim_score = (float(total_sim_score)/num_gens)
avg_sim_aa_edits = (float(total_sim_aa_edits)/num_gens)
sig_pvals = rmath.sig_pvalue_frequency(p_values)
m_o.write("%s,%s,%s,%.2f,%.2f,%.2f,%.2f" % (name,num_edits,num_exp_edits,
    avg_sim_aa_edits,avg_score,avg_sim_score,sig_pvals))
m_o.write("\n")

#print simulation.get_value(simulation.codon_pos_dict)
simulation.write_sim_information(name,s_o)
#simulation.get_codon_percent()
#simulation.get_base_conversion()

# Remove the temporary files when no longer needed
os.remove("tempfile.fa")
os.remove("tempfile.afa")
os.remove("sim_tempfile.fa")
os.remove("sim_tempfile.afa")

# Finally, close the output file again
m_o.close()
