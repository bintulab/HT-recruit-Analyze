###############################################################################
# Josh Tycko
# Associated with HT-recruit paper Tycko et al., 2020

    # Adapted from casTLE CRISPR screen analysis functions by:
    # David Morgens
    # 03/18/2016

###############################################################################
# Imports neccessary modules

'''
Module containing important auxilary screening functions
'''

from __future__ import division
import numpy
import numpy as np
import math
import csv
from collections import defaultdict
import os
import scipy.misc
import scipy.stats as st
import scipy.stats
import sys
import random
import re


###############################################################################
# Function to round off significant digits

def sigDig(x, num=3):
    '''
    Function rounds number in a reasonable manner. Default is to three
    significant digits.
    '''
    # Don't attempt to calculate the log of zero
    if x == 0:
        return 0

    else:
        order_of_magnitude = int(math.floor(math.log10(abs(x))))
        digits = (num - 1) - order_of_magnitude

    return round(x, digits)


# #########################################################################
# # Calculates the GC content of an inputed string

# def getGC(guide):

#     numGC = 0
#     total = len(guide)

#     for nuc in list(guide):
#         if nuc == 'C' or nuc == 'G' or nuc == 'c' or nuc == 'g':
#             numGC += 1

#     return float(numGC) / total

###############################################################################
# Processes time zero count files

def timeZero(zero_files, thresh):

    # Defaults values to count threshold
    zero_unt = defaultdict(lambda: thresh)
    zero_trt = defaultdict(lambda: thresh)

    # If no time zero file provided, returns defaults only
    if not zero_files:
        return zero_unt, zero_trt
    else:
        zero_unt_file, zero_trt_file = zero_files

    # Reads and filters in time zero untreated file
    with open(zero_unt_file, 'rU') as zero_unt_open:

        dialect = csv.Sniffer().sniff(zero_unt_open.read(1024), delimiters='\t ,')
        zero_unt_csv = csv.reader(zero_unt_open, dialect)
        #zero_unt_csv = csv.reader(zero_unt_open, delimiter=',')

        for line in zero_unt_csv:

            if int(line[1]) > thresh:
                zero_unt[line[0]] = int(line[1])

            else:
                zero_unt[line[0]] = thresh

    # Reads and filters in time zero treated file
    with open(zero_trt_file, 'rU') as zero_trt_open:

        dialect = csv.Sniffer().sniff(zero_trt_open.read(1024), delimiters='\t ,')
        zero_trt_csv = csv.reader(zero_trt_open, dialect)
        #zero_trt_csv = csv.reader(zero_trt_open, delimiter=',')

        for line in zero_trt_csv:

            if int(line[1]) > thresh:
                zero_trt[line[0]] = int(line[1])

            else:
                zero_trt[line[0]] = thresh

    return zero_unt, zero_trt


###############################################################################
# Filters count file by a threshold. If counts are below threshold, redefines
# them to equal threshold.  If counts in both samples are below threshold,
# throws them out.

def filterCounts(unt_file, trt_file, thresh, zero_files, exclude=False):
    '''
    Takes untreated and treated count files and filters them according
    to threshold.
    '''

    # Processes time zero files in auxilary function
    zero_unt_raw, zero_trt_raw = timeZero(zero_files, thresh)

    # Stores untreated counts as dictionary of name to count
    untreated_raw = {}
    treated_raw = {}

    with open(unt_file, 'rU') as unt_open:

        dialect = csv.Sniffer().sniff(unt_open.read(1024), delimiters='\t ,')
        unt_open.seek(0)
        unt_csv = csv.reader(unt_open, dialect)
        #unt_csv = csv.reader(unt_open, delimiter=',')

        for line in unt_csv:

            # Skips blank lines
            if not line or not line[0]:
                continue

            # If no exclusion characters, save line
            if not exclude:
                    untreated_raw[line[0]] = int(float(line[1]))

            # If exclusion character if it does not contain substring
            else:
                for ex in exclude:
                    if ex in line[0]:
                        untreated_raw[line[0]] = int(line[1])
                        break

    # Stores treated counts as dictionary of name to count
    with open(trt_file, 'rU') as trt_open:

        dialect = csv.Sniffer().sniff(trt_open.read(1024), delimiters='\t ,')
        trt_open.seek(0)
        trt_csv = csv.reader(trt_open, dialect)
        #trt_csv = csv.reader(trt_open, delimiter=',')

        for line in trt_csv:

            # Skips blank lines
            if not line or not line[0]:
                continue

            # If no exclusion characters, save line
            if not exclude:
                treated_raw[line[0]] = int(float(line[1]))

            # If exclusion character if it does not contain substring
            else:
                for ex in exclude:
                    if ex in line[0]:
                        treated_raw[line[0]] = int(line[1])
                        break         

    # Tracks some filtering statistics
    belowUnt, belowTrt = 0, 0
    removed = 0

    # Stores filtered counts as dictionary of name to count
    treated = {}
    untreated = {}
    zero_unt = {}
    zero_trt = {}

    # Loops over untreated counts, looks for that entry in the the treated
    # counts. Nonpresence indicates zero counts.  If both counts are less
    # than the threshold, the entry is filtered out.  If one is less, then
    # it is assigned to the threshold value. Elsewise saves the values.
    for entry in untreated_raw:

        # Indicator variable of meeting the threshold in each count file
        un = 0
        tr = 0

        # Checks if over threshold in untreated sample
        if untreated_raw[entry] < thresh:
            un = 1
            belowUnt += 1

        # Checks if over threshold in treated sample
        if entry not in treated_raw or treated_raw[entry] < thresh:
            tr = 1
            belowTrt += 1

        # If under in both, don't save the entry
        if un and tr:
            removed += 1
            continue

        # If under threshold in untreated, save as threshold value
        if un:
            untreated[entry] = thresh
            
        else:
            untreated[entry] = untreated_raw[entry]

        # If under threshold in treated, save as threshold value
        if tr:
            treated[entry] = thresh
        else:
            treated[entry] = treated_raw[entry]

        # Looks up time zero counts
        zero_unt[entry] = zero_unt_raw[entry]
        zero_trt[entry] = zero_trt_raw[entry]

    # Loops over treated, looking for entries missed in the untreated counts
    for entry in treated_raw:
        if entry not in untreated_raw:

            # If too small in both, do not save.
            if treated_raw[entry] < thresh:
                removed += 1

            # Else save with untreated value equal to threshold
            else:
                treated[entry] = treated_raw[entry]
                untreated[entry] = thresh
                belowUnt += 1

                # Looks up time zero counts
                zero_unt[entry] = zero_unt_raw[entry]
                zero_trt[entry] = zero_trt_raw[entry]

    # Saves stats and time zero files
    stats = (belowTrt, belowUnt, removed)
    time_zero = (zero_unt, zero_trt)

    return untreated, treated, stats, time_zero
 

###############################################################################
# Function to calculate enrichment values

def enrich(count1, sum1, zero1, sum_zero1,
           count2, sum2, zero2, sum_zero2,
           shift, norm, bkgrd):
    '''
    Function calculates enrichment values
    '''

    # Calculates proportions
    prop1 = float(count1) / sum1
    prop2 = float(count2) / sum2
    prop_zero1 = float(zero1) / sum_zero1
    prop_zero2 = float(zero2) / sum_zero2

    # Calculates ratio and log ratio 
    log_enrich = math.log(prop1 / prop2, 2) - math.log(prop_zero1 / prop_zero2, 2)

    # Normalizes log ratio by shifting around 'zero' and stretching
    # appropriately
    shift_enrich = log_enrich - shift
    norm_enrich = shift_enrich / nor

    return [norm_enrich, count1, count2]


###############################################################################
# Function to calculate enrichments

def enrich_all(untreated, treated, neg_name, split_mark, K, time_zero, back):
    '''
    Auxilary function to calculate enrichment values
    '''

    # Finds total counts in time zero files
    zero_unt, zero_trt = time_zero
    total_zero_unt = sum(zero_unt.values())
    total_zero_trt = sum(zero_trt.values())

    # Finds total counts in count files
    total_unt = sum(untreated.values())
    total_trt = sum(treated.values())

    # Stores enrichments of negative controls
    neg_raw = []
    tar_raw = []
    back_raw = []
    neg_Poff = []

    # Moves over each element, and, if it is a control element, calculates its enrichment
    for entry in untreated:

        # Select only negative controls
        if entry.split(split_mark)[0].startswith(neg_name):

            # Calls enrichment function with a 0 shift
            neg_raw.append(enrich(treated[entry], total_trt,
                                    zero_trt[entry], total_zero_trt,
                                    untreated[entry], total_unt,
                                    zero_unt[entry], total_zero_unt,
                                    0, 1, 0)[0])

            neg_Poff.append(enrich(treated[entry], total_trt,
                                    zero_trt[entry], total_zero_trt,
                                    untreated[entry], total_unt,
                                    zero_unt[entry], total_zero_unt,
                                    0, 1, 0)[3])

        else:
            tar_raw.append(enrich(treated[entry], total_trt,
                                    zero_trt[entry], total_zero_trt,
                                    untreated[entry], total_unt,
                                    zero_unt[entry], total_zero_unt,
                                    0, 1, 0)[0])

    # options for computing shift for log2 enrichments
    if back == 'neg':
        shift = np.median(neg_raw)  # Calculates the shift as a median

    elif back == 'tar':
        shift = np.median(tar_raw)

    elif back == 'all':
        shift = np.median(neg_raw + tar_raw)

    elif back == 'none':
        shift = 0

    else:
        sys.exit('Unrecognized option for background choice: ' + back)

    # Compute background for percent OFF
    bkgrd = np.median(neg_Poff)

    entry_rhos = {}

    # With the shift in hand, calculates the enrichment for each element
    for entry in treated:

        # Note the calculated neg_shift is used now
        entry_rhos[entry] = enrich(treated[entry], total_trt,
                                    zero_trt[entry], total_zero_trt,
                                    untreated[entry], total_unt,
                                    zero_unt[entry], total_zero_unt,
                                    shift, K, bkgrd)

    # Gathers rho values for each gene and separates out negative controls
    gene_rhos = defaultdict(list)
    gene_rhos_int = defaultdict(list)
    gene_ref = defaultdict(list)

    # Stores all negative element rhos and targeting rhos
    neg_rhos = []
    tar_rhos = []

    for entry in entry_rhos:

        # Checks if entry is a negative control
        if entry.split(split_mark)[0].startswith(neg_name):
            neg_rhos.append(entry_rhos[entry])

        else:

            # Gathers rhos of elements targeting each gene
            gene = entry.split(split_mark)[0].upper()
            gene_rhos[gene] += [entry_rhos[entry]]

            # Saves element name and enrichment for output
            entry_split = entry.split(split_mark)
            gene_ref[gene] += [(entry_rhos[entry], entry)]
            tar_rhos.append(entry_rhos[entry])

    return entry_rhos, gene_rhos, neg_rhos, tar_rhos, gene_ref


###############################################################################
# Script to process result records

def retrieveRecord(res_file, current_version):

    name = res_file[: -4]
    rec_file = name + '_record.txt'

    try:
        # Parses record file
        with open(rec_file, 'r') as rec_open:
            rec_csv = csv.reader(rec_open, delimiter='\t')
            script, version = rec_csv.next()

            if version != current_version:
                sys.exit('Error: Version number not current\n'
                            + 'Rerun analysis')

            if script != 'analyzeCounts.py':
                sys.exit('Error: Input is not a result file')

            last_time = rec_csv.next()[1]
            unt_file = rec_csv.next()[1]
            trt_file = rec_csv.next()[1]
            zero_files = rec_csv.next()[1]
            if zero_files:
                zero_files = eval(zero_files)
            file_out = rec_csv.next()[1]
            screen_type = rec_csv.next()[1]
            neg_name = rec_csv.next()[1]
            split_mark = rec_csv.next()[1]
            exclude = rec_csv.next()[1]
            if exclude:
                exclude = eval(exclude)
            thresh = int(rec_csv.next()[1])
            K = float(rec_csv.next()[1])
            back = rec_csv.next()[1]
            I_step = float(rec_csv.next()[1])
            scale = int(rec_csv.next()[1])
            draw_num = int(rec_csv.next()[1])

    except IOError:
        sys.exit('Error: Record of result file not found\n'
                    + 'Change file name or rerun analysis')

    # Saves parameters for processing
    stats = (script, version, last_time)
    files = (unt_file, trt_file, zero_files, file_out)
    info = (screen_type, neg_name, split_mark, exclude)
    param = (thresh, K, back, I_step, scale, draw_num)

    return stats, files, info, param

###############################################################################
