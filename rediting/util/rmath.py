#!/usr/bin/env python

"""This module includes functions to calculate various metrics. This
includes simple calculations such as mean, but also specific calcs for
things like pearson correlation coefficient and polyT"""

import math
import collections
import random
from bisect import bisect
import numpy as np

from rediting.classes import matrices
import sequence

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
        if sequence.polyTpercent(e, percent):
            # If at least one window is True, than evaluates True
            return True
        else:
            pass
    # This only returns False if ALL windows are False
    return False


def calculate_entropy(aa_list):
    """Calculates entropy for a list of values"""
    # Calculate aa entropy without gaps
    new_aa_list = [x for x in aa_list if x != '-']
    num_gaps = len(aa_list) - len(new_aa_list)
    # Determine how much of the total position is gaps
    gap_proportion = float(num_gaps)/len(aa_list)
    num_aa = len(new_aa_list)
    aa_count = collections.Counter(new_aa_list)
    total_entropy = 0
    for k,v in aa_count.items():
        aa_proportion = (float(v)/num_aa)
        # Calculate normalized entropy
        entropy = aa_proportion * math.log(aa_proportion) * (1/math.log(20))
        total_entropy += entropy
    total_entropy = (total_entropy * -1)
    # End result takes into account both gaps and entropy
    result = (1 - gap_proportion) * (1 - total_entropy)
    return result


# Might not need this after all
def calculate_score_diff(aa_edit_list):
    """calculates the average score difference based on
    a provided list of aa editing characteristics"""
    num_edits = float(len(aa_edit_list))
    diff_score = 0
    for P,RA,GA,MA,IB,IA,SB,SA,D in aa_edit_list:
        diff_score += int(D)
    try:
        avg_diff = diff_score/num_edits
    except(ZeroDivisionError):
        avg_diff = 0
    return avg_diff


def calculate_average_score_diff(rg,rm,aa_list):
    """Calculates average edit score difference between
    genomic/transcript amino acid and reference sequences"""
    total_score = 0
    for ref_aa in aa_list:
        scr_before = matrices.Blosum62(rg,ref_aa).sub_score()
        scr_after = matrices.Blosum62(rm,ref_aa).sub_score()
        total_score += (scr_after - scr_before)
    return (float(total_score)/len(aa_list))


def weighted_choice(choices):
    """Allows random choices but with weighted values. Takes a list of
    2-length tuples as input. Code based on:

    http://stackoverflow.com/questions/3679694/a-weighted-version-of-random-choice
    """
    values, weights = zip(*choices)
    total = 0
    cumulative_weight = []
    for w in weights:
        total += w
        cumulative_weight.append(total)
    rand_val = random.random() * total
    i = bisect(cumulative_weight, rand_val)
    return values[i]


def equal_vars(array1,array2):
    """Checks whether the variance is equal"""
    if np.var(array1) == np.var(array2):
        return True
    else:
        return False

def sig_pvalue_frequency(p_values,threshold=0.05):
    """Calculates the proportion of p_values in a list
    that are below a given significance threshold"""
    sig_pvals = 0
    for p_val in p_values:
        if p_val < threshold:
            sig_pvals += 1
    try:
        return (float(sig_pvals)/len(p_values))
    except(ZeroDivisionError):
        return 0.0
