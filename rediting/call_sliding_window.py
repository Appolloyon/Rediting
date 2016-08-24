#!/usr/bin/env python

import subprocess
import argparse

parser = argparse.ArgumentParser(
    description = "Calls the sliding window program for one or more files",
    epilog = """This program allows for batch calls of the sliding window
    program. If files have routine structure and are named in a consistent
    manner including information such as organism and gene name then a simple
    wrapper script such as this can easily automate largescale analysis""")
parser.add_argument('infiles', nargs='+', help='list of infiles')
parser.add_argument('-d', '--delimiter', help='delimiter to separate file info')
parser.add_argument('-g', '--genefield', help='field for gene name from 1 to X')
parser.add_argument('-gen', '--genomic', help='unique string present in genomic sequenc headers')
parser.add_argument('-r', '--RNA', help='unique string present in RNA sequence headers')
parser.add_argument('-out', '--outfile', help='name for master outfile')
args = parser.parse_args()

# Determines how to split filenames
delimiter = args.delimiter
# Given delimeter, which filename field has the gene name
genefield = int(args.genefield) - 1
gen = args.genomic
rna = args.RNA
out = args.outfile

for infile in args.infiles:
    # Gets the basename of the file
    short_in = infile.rsplit('.',1)[0]
    # Next two lines get the gene name
    filename_list = short_in.split(delimiter)
    gene = filename_list[genefield]
    # For each file, call the program
    subprocess.call(["/Users/cklinger/git/Rediting/rediting/sliding_window.py",
        "-in", infile, "-out", out, "-n", short_in,
        "-g", gene, "-r", rna, "-gen", gen, "-p"])#, "-b", "-o"])
