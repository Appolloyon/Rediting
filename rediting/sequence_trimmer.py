#!/usr/bin/env python

import re
import argparse

from classes import classes
from util import files,strings,sequence

parser = argparse.ArgumentParser(
    description = """Trims aligned genomic/RNA sequences and a reference
        sequence to the minimum possible length""",
    epilog = """This program assumes that the genomic and RNA sequence for a
    given gene are provided in the same file in an aligned FASTA format, such
    as that output by MAFFT or MUSCLE. In addition, a reference sequence in
    the same aligned file will be used for comparison. This program will step
    through the aligned sequences and check for regions present in both the
    reference and genomic sequence, and remove as many of these as possible
    while maintaining the reading frame of each sequence. This will allow for
    calculation with the sliding window program without the potential artifact
    that long stretches of gaps in either sequence may introduce.""")
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequence header')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence header')
parser.add_argument('infiles', nargs='+', help='list of infiles')
args = parser.parse_args()

for infile in args.infiles:
    basename = infile.rsplit('.',1)[0] #everything before final period
    outfile = basename + "_trimmed.fa"

    seqdict = {}
    files.build_seqdict(infile,seqdict)

    rna_string = str(args.RNA)
    gen_string = str(args.genomic)
    # Sequences must be in upper case
    for k in seqdict.keys():
        if re.search(rna_string,k):
            rna_header = k
            rna_seq = seqdict.get(k).upper()
        elif re.search(gen_string,k):
            gen_header = k
            gen_seq = seqdict.get(k).upper()
        else:
            ref_header = k
            ref_seq = seqdict.get(k).upper()

    san_gen_seq = strings.sanitize(gen_seq)
    san_ref_seq = strings.sanitize(ref_seq)
    ref_pair = classes.RefPair(san_ref_seq,san_gen_seq,"name")

    gen_start = 'NA' # can't use False, as zero index also evaluates
    ref_start = 'NA'
    gen_list = []
    ref_list = []
    for i, (rg,rf) in enumerate(zip(gen_seq,ref_seq)):
        #print "%d, rg is %s, rf is %s, start is %s" % (i,rg,rf,ref_start)
        if rf == '-' and rg == '-':
            pass
        # Identify inserts in reference sequence
        elif rf != '-' and rg == '-' and ref_start == 'NA':
            if ref_pair.index_rposition() == 1:
                ref_start = i
            else:
                pass
            ref_pair.incr_all_ref()
        elif rf != '-' and rg == '-' and ref_start != 'NA':
            if i == (len(ref_seq) - 1):
                ref_end = sequence.close_gap(i,ref_pair.index_rposition())
                ref_list.append([ref_start,ref_end])
                ref_pair.incr_all_ref()
                ref_start = 'NA'
            else:
                ref_pair.incr_all_ref()
        elif rf != '-' and rg != '-' and ref_start != 'NA':
            ref_end = sequence.close_gap(i,ref_pair.index_rposition())
            ref_list.append([ref_start,ref_end])
            ref_pair.incr_all_ref()
            ref_start = 'NA'
        # Repeat for genomic sequence
        elif rg != '-' and rf == '-' and gen_start == 'NA':
            if ref_pair.index_gposition() == 1:
                gen_start = i
            else:
                pass
            ref_pair.incr_all_gen()
        elif rg != '-' and rf == '-' and gen_start != 'NA':
            if i == (len(gen_seq) - 1):
                gen_end = sequence.close_gap(i,ref_pair.index_gposition())
                gen_list.append([gen_start,gen_end])
                ref_pair.incr_all_gen()
                gen_start = 'NA'
            else:
                ref_pair.incr_all_gen()
        elif rg != '-' and rf != '-' and gen_start != 'NA':
            gen_end = sequence.close_gap(i,ref_pair.index_gposition())
            gen_list.append([gen_start,gen_end])
            ref_pair.incr_all_gen()
            gen_start = 'NA'


    ref_indices = []
    for (start,stop) in ref_list:
        if sequence.check_indices(start,stop):
            sequence.expand_indices(start,stop,ref_indices)

    gen_indices = []
    for (start,stop) in gen_list:
        if sequence.check_indices(start,stop):
            sequence.expand_indices(start,stop,ref_indices)

    new_ref_seq = sequence.trim_sequence(ref_seq,ref_indices)
    new_gen_seq = sequence.trim_sequence(gen_seq,gen_indices)
    new_rna_seq = sequence.trim_sequence(rna_seq,gen_indices)

    with open(outfile,'w') as o:
        o.write(">" + ref_header + "\n")
        for chunk in strings.split_input(new_ref_seq,60):
            o.write(chunk + "\n")
        o.write(">" + gen_header + "\n")
        for chunk in strings.split_input(new_gen_seq,60):
            o.write(chunk + "\n")
        o.write(">" + rna_header + "\n")
        for chunk in strings.split_input(new_rna_seq,60):
            o.write(chunk + "\n")
