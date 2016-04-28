#!/usr/bin/env python

import argparse

from util import sequence

parser = argparse.ArgumentParser(
    description = """Generates random nucleotide sequences for comparison""",
    epilog = """This program is intended to generate a series of nucleotide sequence
    pairs, wherein the starting sequence is always the same for each program run but
    each iteration produces a randomly mutated version.""")
parser.add_argument('-l', '--length', help='length of sequence')
parser.add_argument('-n', '--number', help='number of mutations')
parser.add_argument('-g', '--generations', help='number of pairs to generate')
args = parser.parse_args()

length = int(args.length)
num_muts = int(args.number)
num_gens = int(args.generations)
bases = 'AGTC'

starting_sequence = sequence.generate_start_sequence(length,bases)

i = 0
while i < num_gens:
    out_name = "seq_pair_" + str(i+1) + ".fa"
    with open(out_name,'w') as o:
        new_seq = sequence.mutate_sequence(starting_sequence,num_muts,bases)
        o.write(">original_seq" + "\n" + "".join(starting_sequence) + "\n")
        o.write(">new_seq" + "\n" + "".join(new_seq) + "\n")
    i += 1
