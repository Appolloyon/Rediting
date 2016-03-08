"""."""

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
