#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser(
    description = """Generates random mutated versions of a genomic sequence,
        based on known codon position and base conversion biases"""
    epilog = """This program takes a file with a single sequence and parameters
    describing the known biases in editing for the relevant organism and
    simulates the effect of editing on the sequence for a number of specified
    edits and any number of independent times. Output is provided as a series
    of files, each named according to a specified convention.""")
parser.add_argument('-in', '--infile', help='file with starting sequence')
parser.add_argument('-n', '--number', help='number of edits to simulate')
parser.add_argument('-g', '--generations', help='number of simulations')
parser.add_argument('-c', '--codon', help='file with codon bias information')
parser.add_argument('-x', '--exchange', help='file with base convserion information')
args = parser.parse_args()




"""
Pseudo-code for program execution

Supply information for execution:
    -sequence to use - this can simply be the string
    -number of edits - this is just a value
    -distribution of codon edits - should probably read it from a file
    -distribution of base edits - should definitely read it from a file

Partition the sequence into codon positions
    -make three separate lists of positions, indexed by value

Mutate the sequence, which means for each edit:
    -choose whether it will be in the 1st, 2nd, or 3rd codon position
        -use the weighted choice function
    -get a random position from the corresponding list
    -check what the nucleotide is at that position
    -mutate it with a weighted choice

Write the new sequence to an output file
"""
