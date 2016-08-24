#!/usr/bin/env python

"""This module contains functions pertaining to dealing with reading from,
and writing to, files as well as organizing the resulting data from read
files into usable data structures."""


def nonblank_lines(f):
    """read file lines, but skip blank lines"""
    for l in f:
        line = l.strip('\n')
        if line: # Equivalent to if 'line not blank'
            yield line


def build_seqdict(infile, seqdict):
    """Builds a dictionary of header/sequence key value pairs.
    NOTE: assumes FASTA formatted input file!"""
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


def get_variable_file_handle(infile,mode,delimeter=None,write_list=None):
    """Returns an open file object. If write_list is specified,
    the individual elements of the list are written to the file
    using delimeter as a separator"""
    file_obj = open(infile,mode)
    if write_list is not None and delimeter is None:
        # If not specified, but needed, delimeter is a comma
        delimeter = ','
    # Typically, use this to write header lines
    if write_list is not None:
        for elem in write_list:
            file_obj.write(elem + delimeter)
        file_obj.write("\n" * 2)
    return file_obj


def write_transition_dict(gene,tdict,file_obj):
    """Takes a dictionary and file object and writes the
    relevant information to the file"""
    file_obj.write("%s,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d,%d" % (gene,
        tdict.get('a_t'),tdict.get('a_g'),tdict.get('a_c'),
        tdict.get('t_a'),tdict.get('t_g'),tdict.get('t_c'),
        tdict.get('g_a'),tdict.get('g_t'),tdict.get('g_c'),
        tdict.get('c_a'),tdict.get('c_t'),tdict.get('c_g')))
    file_obj.write("\n")


# Potentially obsolete
def parse_codon_bias(infile, codon_weights):
    """Builds a list of codon position bias weights"""
    with open(infile,'U') as f:
        for line in nonblank_lines(f):
            llist = line.strip('\n').split(',')
            for i,w in enumerate(llist):
                codon_weights.append(((i+1),float(w)))
    return codon_weights


# Potentially obsolete
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


def parse_codon_base_bias(infile,weights):
    """Function definition to go here"""
    with open(infile,'U') as f:
        for line in nonblank_lines(f):
            llist = line.strip('\n').split(',')
            codon_pos = llist[0]
            start_base = llist[1]
            end_base = llist[2]
            freq = int(llist[3])
            weights.append(((codon_pos,start_base,end_base),freq))
    return weights
