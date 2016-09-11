#!/usr/bin/env python

import subprocess
import argparse

parser = argparse.ArgumentParser(
    description = "Calls the entropy program for one or more files",
    epilog = """This program allows for batch calls of the entropy calc
    program. If files have routine structure and are named in a consistent
    manner including information such as organism and gene name then a simple
    wrapper script such as this can easily automate largescale analysis""")
parser.add_argument('infiles', nargs='+', help='list of infiles')
parser.add_argument('-d', '--delimiter', help='delimiter to separate file info')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequenc headers')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence headers')
args = parser.parse_args()

# Deteremines how to split filenames
delimiter = args.delimiter
# Given delimeter, which filename field has the gene name
gen = args.genomic
rna = args.RNA

for infile in args.infiles:
    # Gets the basename of the file
    short_in = infile.rsplit('.',1)[0]
    # For each file, call the program
    subprocess.call(["/Users/cklinger/git/Rediting/rediting/calculate_entropy.py",
            "-in", infile, "-n", short_in, "-r", rna, "-gen", gen])
