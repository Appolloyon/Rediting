"""."""

import math


def calc_mean(values):
    """Calculates mean of a set of values"""
    sum = 0.0
    for value in values:
        value = float(value)
        sum += value
    mean = sum/(float(len(values)))
    return mean


def calc_pearson(xvalues, yvalues, xmean, ymean):
    """Calculates a pearson correlation value"""
    N = len(xvalues)
    num = 0.0
    xdenom = 0.0
    ydenom = 0.0
    for i in range(N):
        num += ((xvalues[i] - xmean) * (yvalues[i] - ymean))
        xdenom += ((xvalues[i] - xmean)**2)
        ydenom += ((yvalues[i] - ymean)**2)
    denom = (math.sqrt(xdenom)) * (math.sqrt(ydenom))
    try:
        return num/denom
    except:
        return 0.0


def calc_tvalue(PC, N):
    """Calculates a t value for a given pearson coefficient"""
    return abs((PC * math.sqrt(N-2))/(math.sqrt(1-(PC**2))))


def ispolyTpercent(plist, percent):
    """check list elements for at least one polyT stretch"""
    for e in plist:
        if polyTpercent(e, percent):
            return True
        else:
            pass
    return False
