#!/usr/bin/env python

import re
import argparse
#from indelclass import IndelPair
from functions import gulp,nonblank_lines,sanitize

parser = argparse.ArgumentParser(
    description = "Calculates position of indels against a ref sequence",
    epilog = """This program assumes that the reference and query sequence for
    a given gene are provided in the same format in aligned FASTA format (.afa).
    It also assumes that the header line for the references sequence will have the
    string 'reference' somewhere. Output is given as a single csv file per gene""")
parser.add_argument('infiles', nargs='+', help='list of aligned files')
args = parser.parse_args()


for infile in args.infiles:
    name = infile.rstrip('.afa')
    with open(infile,'U') as i:
        seqdict = {}
        for line in nonblank_lines(f):
            line = line.strip('\n')
            if line.startswith(">"):
                line = line.lstrip(">")
                ID = line
                seqdict[ID] = ''
            else:
                seqdict[ID] += line

    for k in seqdict.keys():
        if re.search('reference',k):
            rseq = seqdict.get(k)
        else:
            qseq = seqdict.get(k)

    san_rseq = sanitize(rseq)
    san_qseq = sanitize(qseq)

    #seq_pair = IndelPair(san_rseq, san_qseq, name)
    ref_counter = 1

# come back to this later and change the definition of start and end of alignment
    i = 0
    while not (rseq[i] != '-' and qseq[i] != '-'):  # start when both residues exist
        if rseq[i] != '-':  # we only care about the refseq index
            ref_counter += 1
        i += 1
    j = 0
    while not (rseq[(len(rseq)-1)-j] != '-' and qseq[(len(qseq)-1)-j] != '-'):  # same end criterion
        j += 1

    new_rseq = rseq[i:(len(rseq)-j)]
    new_qseq = qseq[i:(len(qseq)-j)]

    insert_dict = {}
    deletion_dict = {}

    for i,(raa,qaa) in enumerate(zip(new_rseq,new_qseq)):
        if raa == '-' and qaa != '-':  # insertion in query
            rpos = ref_counter
            if rpos in insert_dict.keys():
                insert_dict[rpos] += 1
            else:
                insert_dict[rpos] = 1
            ref_counter += 1
        elif raa != '-' and qaa == '-': # deletion in query
            rpos = ref_counter
            deletion_dict[rpos] = 1
        else: # residue in both
            ref_counter += 1
