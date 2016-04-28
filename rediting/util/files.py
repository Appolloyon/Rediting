#!/usr/bin/env python

"""This module contains functions pertaining to dealing with reading from
files and organizing the resulting data into data structures."""


def nonblank_lines(f):
    """read file lines, but skip blank lines"""
    for l in f:
        line = l.strip('\n')
        if line: # Equivalent to if 'line not blank'
            yield line


def build_seqdict(infile, seqdict):
    """Builds a dictionary of header/sequence key value pairs"""
    with open(infile,'U') as f:
        for line in nonblank_lines(f):
            line = line.strip('\n')
            # These are the header lines
            if line.startswith(">"):
                line = line.strip(">")
                ID = line
                seqdict[ID] = ''
            # If not a header line, it must be part of the sequence
            # for the previous header, so add it to the dict value
            else:
                seqdict[ID] += line
    return seqdict


def parse_codon_bias(infile, codon_weights):
    """Builds a list of codon position bias weights"""
    with open(infile,'U') as f:
        for line in nonblank_lines(f):
            llist = line.strip('\n').split(',')
            for i,w in enumerate(llist):
                codon_weights.append(((i+1),float(w)))
    return codon_weights


def parse_base_bias(infile, base_dict, base_weights):
    """Builds a dict of lists for base conversion, and also
    determines the relative frequency of each changing base"""
    # Keep track of all total bases as well
    total_freq_dict = {'A':0,'G':0,'T':0,'C':0}
    with open(infile,'U') as f:
        for line in nonblank_lines(f):
            llist = line.strip('\n').split(',')
            start_base = llist[0]
            end_base = llist[1]
            freq = float(llist[2])
            # Update for each base
            if start_base not in base_dict.keys():
                base_dict[start_base] = [(end_base,freq)]
            else:
                base_dict[start_base].append((end_base,freq))
            # Update for total bases
            total_freq_dict[start_base] += freq
    for k,v in total_freq_dict.items():
        base_weights.append((k,v))
    return base_dict,base_weights

