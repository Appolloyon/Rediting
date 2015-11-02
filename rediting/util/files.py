"""."""


def nonblank_lines(f):
    """skip blank lines"""
    for l in f:
        line = l.strip('\n')
        if line:
            yield line


def build_seqdict(infile, seqdict):
    """Builds a dictionary of header/sequence key value pairs"""
    with open(infile,'U') as f:
        for line in nonblank_lines(f):
            line = line.strip('\n')
            if line.startswith(">"):
                line = line.strip(">")
                ID = line
                seqdict[ID] = ''
            else:
                seqdict[ID] += line
    return seqdict
