#!/usr/bin/env python

"""This module contains string utility functions."""


def gulp(string, start, gulp_size):
    """get substrings of a string"""
    return string[start:start+gulp_size]

def sanitize(seq):
    """remove gap characters"""
    nseq = ''
    for char in seq:
        if char == '-':
            pass
        else:
            nseq += char
    return nseq

def split_input(string, chunk_size):
    """split a string into multiple substrings"""
    num_chunks = len(string)/chunk_size
    if (len(string) % chunk_size != 0):
        # This last chunk is not the same size, but string
        # slicing takes care of the details for us
        num_chunks += 1
    output = []
    for i in range(0, num_chunks):
        output.append(string[chunk_size*i:chunk_size*(i+1)])
    return output
