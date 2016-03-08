"""."""

"""This module includes functions to calculate various metrics. This
includes simple calculations such as mean, but also specific calcs for
things like pearson correlation coefficient and polyT"""

import math

from sequence import polyTpercent

def calc_mean(values):
    """Calculates the mean of a set of values"""
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
        # The numerator corresponds to the difference between
        # each value and its corresponding mean
        num += ((xvalues[i] - xmean) * (yvalues[i] - ymean))
        # Each part of the denominator is squared
        xdenom += ((xvalues[i] - xmean)**2)
        ydenom += ((yvalues[i] - ymean)**2)
    denom = (math.sqrt(xdenom)) * (math.sqrt(ydenom))
    try:
        return num/denom
    # If there is nothing to calculate, then the denominator
    # could theoretically be zero, hence return 0.0 instead
    except:
        return 0.0


def calc_tvalue(PC, N):
    """Calculates a t value for a given pearson coefficient"""
    # Simple calculation based on the pearson value
    return abs((PC * math.sqrt(N-2))/(math.sqrt(1-(PC**2))))


def ispolyTpercent(plist, percent):
    """check list elements for at least one polyT stretch"""
    # Go over each window in a list
    for e in plist:
        if polyTpercent(e, percent):
            # If at least one window is True, than evaluates True
            return True
        else:
            pass
    # This only returns False if ALL windows are False
    return False
