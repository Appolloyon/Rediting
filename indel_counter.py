#!/usr/bin/env python

import re
import argparse

parser = argparse.ArgumentParser(
    description = "Calculates position of indels against a ref sequence",
    epilog = """This program assumes that the reference and query sequence for
    a given gene are provided in the same format in aligned FASTA format.
    It also assumes that the header line for the references sequence will have the
    string 'reference' somewhere. Output is given as a single csv file per gene""")
parser.add_argument('infiles', nargs='+', help='list of aligned files')
args = parser.parse_args()


def nonblank_lines(f):
    """skip blank lines"""
    for l in f:
        line = l.strip('\n')
        if line:
            yield line


for infile in args.infiles:
    name = infile.split('.')[0]
    with open(infile,'U') as i:
        seqdict = {}
        for line in nonblank_lines(i):
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

    ref_counter = 1

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

    insert_list = []
    insert_dict = {}
    deletion_list = []
    deletion_dict = {}

    prev_del = 0 # initialize a counter for the previous deletion

    for i,(raa,qaa) in enumerate(zip(new_rseq,new_qseq)):
        if raa == '-' and qaa != '-':  # insertion in query
            rpos = ref_counter
            if rpos in insert_dict.keys():
                insert_dict[rpos] += 1
            else:
                insert_dict[rpos] = 1
                insert_list.append(rpos)
        elif raa != '-' and qaa == '-': # deletion in query
            rpos = ref_counter
            if abs(rpos - prev_del) == 1:
                key = deletion_list[(len(deletion_list)-1)]
                deletion_dict[key] += 1
                prev_del = rpos
            else:
                deletion_dict[rpos] = 1
                deletion_list.append(rpos)
                prev_del = rpos
            ref_counter += 1
        else: # residue in both
            ref_counter += 1

    out_insert = name + "_inserts.csv"
    out_deletion = name + "_deletions.csv"

    with open(out_insert,'w') as o1, open(out_deletion,'w') as o2:
        for k in insert_list:
            v = insert_dict[k]
            o1.write("%s,%s,%s" % (name,k,v))
            o1.write('\n')
        for k in deletion_list:
            v = deletion_dict[k]
            o2.write("%s,%s,%s" % (name,k,v))
            o2.write('\n')
