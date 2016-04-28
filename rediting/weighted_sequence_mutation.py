#!/usr/bin/env python

import re
import argparse

from classes import classes
from util import files, sequence

parser = argparse.ArgumentParser(
    description = """Generates random mutated versions of a genomic sequence,
        based on known codon position and base conversion biases""",
    epilog = """This program takes a file with a single sequence and parameters
    describing the known biases in editing for the relevant organism and
    simulates the effect of editing on the sequence for a number of specified
    edits and any number of independent times. Output is provided as a series
    of files, each named according to a specified convention.""")
parser.add_argument('-in', '--infile', help='file with starting sequence')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequence headers')
parser.add_argument('-n', '--number', help='number of edits to simulate')
parser.add_argument('-g', '--generations', help='number of simulations')
parser.add_argument('-w', '--weights', help='file with weight values')
#parser.add_argument('-c', '--codon', help='file with codon bias information')
#parser.add_argument('-x', '--exchange', help='file with base convserion information')
args = parser.parse_args()

# Initialize important variables
num_muts = int(args.number)
num_gens = int(args.generations)
gen_string = str(args.genomic)

# Get the genomic sequence
seqdict = {}
files.build_seqdict(args.infile,seqdict)
for k,v in seqdict.items():
    if re.search(gen_string,k):
        gen_seq = seqdict.get(k).upper()
# Transform sequence into a list for future use - need mutable object
gen_list = []
for b in gen_seq:
    gen_list.append(b)

# Create indices for each codon position in sequence
codon_positions,codon_bases = sequence.get_codon_position_indices(gen_seq)
#print codon_bases
# Parse input files to gather frequency information
#codon_weights = []
#files.parse_codon_bias(args.codon,codon_weights)
#base_dict = {}
#base_weights = []
# We need information on relative frequency of change for each base
# and base-specific conversion rates
#files.parse_base_bias(args.exchange,base_dict,base_weights)
weights = []
files.parse_codon_base_bias(args.weights,weights)


# Mutate sequence and write to output files
#i = 0
#while i < num_gens:
#    out_name = "placeholder_" + str(i+1) + ".fa"
#    with open(out_name,'w') as o:
#        new_seq = sequence.weighted_mutation(gen_list,num_muts,
#            codon_positions,codon_weights,base_dict,base_weights)
#        o.write(">new_seq_" + str(i+1) + "\n" + "".join(new_seq) + "\n")
#    i += 1

#sequence.get_overall_weights(codon_weights,codon_bases,base_weights)

simulation = classes.Simulation()
i = 0
while i < num_gens:
    new_seq = sequence.weighted_mutation_new(gen_list,num_muts,
        codon_positions,weights,simulation)
    i += 1

simulation.get_codon_percent()
print
simulation.get_base_conversion()
