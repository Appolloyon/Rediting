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

# Unlike other programs in this package, this one is written to be used
# without a wrapper script, but could be easily adapted to do so
for infile in args.infiles:
    # Gets the basename for the file
    basename = infile.rsplit('.',1)[0]
    # We actually provide aligned and sequence-only
    # versions of the output
    out_align = basename + "_trimmed.afa"
    out_seq = basename + "_trimmed.fa"

    # Load sequence data into a data structure for internal use
    seqdict = {}
    files.build_seqdict(infile,seqdict)

    rna_string = str(args.RNA)
    gen_string = str(args.genomic)
    # Sequences must be in upper case
    for k in seqdict.keys():
        if re.search(rna_string,k):
            # Since we are writing these data back out again
            # we want to keep track of sequence headers
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
    # We steal the RefPair class, but we do not care about the name for
    # writing to output, use "name" as a placeholder
    ref_pair = classes.RefPair(san_ref_seq,san_gen_seq,"name")

    gen_start = 'NA' # Can't use False, as zero index also evaluates
    ref_start = 'NA'
    gen_list = []
    ref_list = []
    for i, (rg,rf) in enumerate(zip(gen_seq,ref_seq)):
        # A gap in both genomic and reference sequences is unlikely,
        # but we should account for it just in case
        if rf == '-' and rg == '-':
            pass

        # Identify inserts in reference sequence
        # We found an indel in the ref relative to the genomic
        elif rf != '-' and rg == '-' and ref_start == 'NA':
            # In order to maintain the reading frame, we actually need
            # to check the codon position
            if ref_pair.index_rposition() == 1:
                # If our gap starts in RF +1 then we can continue
                ref_start = i
            else:
                # If not, keep checking
                pass
            ref_pair.incr_all_ref()
        # We continue an indel in the reference sequence
        elif rf != '-' and rg == '-' and ref_start != 'NA':
            # If the indel is at the end, then we need to close it
            if i == (len(ref_seq) - 1):
                ref_end = sequence.close_gap(i,ref_pair.index_rposition())
                ref_list.append([ref_start,ref_end])
                ref_pair.incr_all_ref()
                ref_start = 'NA'
            else:
                # If not at the end, keep going
                ref_pair.incr_all_ref()
        # We started an indel and now we find a residue in the genomic sequence
        elif rf != '-' and rg != '-' and ref_start != 'NA':
            # Simply close the gap
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


    # Check over all of the reference indices
    ref_indices = []
    for (start,stop) in ref_list:
        # Check to make sure that indices are actually realistic
        # They will not be added if they are less than 3 bases, if
        # they are not multiples of three, or if they are invalid
        if sequence.check_indices(start,stop):
            # Add all of the individual index values to another list
            sequence.expand_indices(start,stop,ref_indices)

    # Check over all of the genomic indices
    gen_indices = []
    for (start,stop) in gen_list:
        if sequence.check_indices(start,stop):
            sequence.expand_indices(start,stop,gen_indices)

    # Technically, we need to remove all non-redundant ones from each
    indices = list(set(ref_indices) | set(gen_indices))

    # Finally, trim the sequences
    new_ref_seq = sequence.trim_sequence(ref_seq,indices)
    new_gen_seq = sequence.trim_sequence(gen_seq,indices)
    # It is assumed that the rna sequence is similar to the genomic
    # sequence, and is trimmed using the indices determined by the
    # other two sequences
    new_rna_seq = sequence.trim_sequence(rna_seq,indices)

    # Write aligned sequences
    with open(out_align,'w') as o:
        o.write(">" + ref_header + "\n")
        for chunk in strings.split_input(new_ref_seq,60):
            o.write(chunk + "\n")
        o.write(">" + gen_header + "\n")
        for chunk in strings.split_input(new_gen_seq,60):
            o.write(chunk + "\n")
        o.write(">" + rna_header + "\n")
        for chunk in strings.split_input(new_rna_seq,60):
            o.write(chunk + "\n")

    # Write sequences without gaps
    with open(out_seq,'w') as o2:
        o2.write(">" + ref_header + "\n")
        for chunk in strings.split_input(strings.sanitize(new_ref_seq),60):
            o2.write(chunk + "\n")
        o2.write(">" + gen_header + "\n")
        for chunk in strings.split_input(strings.sanitize(new_gen_seq),60):
            o2.write(chunk + "\n")
        o2.write(">" + rna_header + "\n")
        for chunk in strings.split_input(strings.sanitize(new_rna_seq),60):
            o2.write(chunk + "\n")
