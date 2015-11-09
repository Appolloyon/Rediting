"""String utility functions."""


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
