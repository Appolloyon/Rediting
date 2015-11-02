"""."""


def gulp(string, start, gulp_size):
    """get substrings of a string"""
    gulpstr = ''
    chars = string[start:start+gulp_size]
    for char in chars:
        gulpstr += char
    return gulpstr


def sanitize(seq):
    """remove gap characters"""
    nseq = ''
    for char in seq:
        if char == '-':
            pass
        else:
            nseq += char
    return nseq
