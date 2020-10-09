###############################################################################
# Josh Tycko
# 03/25/2019
# Adapted from casTLE (by David Morgens)
###############################################################################
# Import neccessary modules

'''
Function for comparing count files.
'''

from __future__ import division
import csv
import time
import argparse
from screenFun import *
import sys
import numpy as np


###############################################################################    
# Version number

current_version = '1.0'

###############################################################################    

# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Compares count files using casTLE')

# Non-optional arguments:
parser.add_argument('ON_file', help='File for ON counts', type=str)

parser.add_argument('OFF_file', help='File for OFF counts', type=str)

parser.add_argument('name', help='Name for output files', type=str)

# Options for element IDs
parser.add_argument('-n', '--neg_name',
                help='Symbol used to denote negative controls',
                type=str, default='0')

parser.add_argument('-s', '--split', dest='split_mark', help='Delimiter for element name',
                type=str, default='_')

parser.add_argument('-x', '--exclude', type=str, help='Excludes substrings', nargs='+')

# Optional arguments for how analysis is performed:
parser.add_argument('-t', '--threshhold', dest='thresh',
                help='Read cutoff for small count numbers', type=int, default=5)

parser.add_argument('-k', '--strength', dest='K', help='Normalizing constant',
                type=float, default=1.0)

parser.add_argument('-b', '--back', help='Background population for noise estimation',
                default='neg', choices=['all', 'neg', 'tar', 'none'])

parser.add_argument('-z', '--zero_files', help='Time zero count files',
                nargs=2, type=str, default='')

# Options for how computation is performed
parser.add_argument('-c', '--scale', type=int, default=3,
                help='Scale of calculations; default is 3')

parser.add_argument('-I', '--I_step', type=float, default=0.1,
                help='Step size in grid search; default is 0.1')

parser.add_argument('-p', '--proccessors', dest='nums',
                help='Number of proccessors to use; default is 20', type=int,
                default=20)

# Options for overriding behavior
parser.add_argument('-r', '--reference',
                help='Location of reference files; default is GeneRef', type=str,
                default='GenRef')

parser.add_argument('-of', '--override_file', action='store_true',
                help='Overrides restriction of output to Results folder')

parser.add_argument('-m', '--mouse', action='store_true',
                help='Uses mouse gene information')

parser.add_argument('-ro', '--record', action='store_false',
                help='Allows script to run without record of count files')

parser.add_argument('-a', '--add', type=str, default='')

# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Processes and checks input arguments

# Creates output file name
if args.override_file:
    file_out = args.name
else:
    file_out = os.path.join('Data', args.name)

# Checks if can write to output
try:
    with open(file_out + '_record.txt', 'w') as out_open:
        pass

    os.remove(file_out + '_record.txt')

except:
    sys.exit('Cannot write to output file:\n' + file_out + '\n'
                + 'Use -of or --override_file to change')


###############################################################################
# Finds record files

if args.record:

    print('Retrieving records')

    # Retrieves record files for count files
    ON_rec_name = args.ON_file[: -11]
    ON_rec_file = ON_rec_name + '_record.txt'

    OFF_rec_name = args.OFF_file[: -11]
    OFF_rec_file = OFF_rec_name + '_record.txt'

    try:
        # Parses record files
        with open(ON_rec_file, 'r') as rec_open:
            rec_csv = csv.reader(rec_open, delimiter='\t')
            ON_version = next(rec_csv)[1]
            ON_time1 = next(rec_csv)[1]
            ON_seq_file = next(rec_csv)[1]
            ON_seq_add_file = next(rec_csv)[1]
            ON_out = next(rec_csv)[1]
            ON_screen_type = next(rec_csv)[1]

        with open(OFF_rec_file, 'r') as rec_open:
            rec_csv = csv.reader(rec_open, delimiter='\t')
            OFF_version = next(rec_csv)[1]
            OFF_time1 = next(rec_csv)[1]
            OFF_seq_file = next(rec_csv)[1]
            OFF_seq_add_file = next(rec_csv)[1]
            OFF_out = next(rec_csv)[1]
            OFF_screen_type = next(rec_csv)[1]

    except IOError:

        print(ON_rec_file)
        print(OFF_rec_file)

        sys.exit('Record of count file not found\n'
                + 'Change file name or rerun makeCounts.py')

    # Checks comparison makes sense
    if ON_screen_type != OFF_screen_type:
        sys.exit('Screen types do not match.')

else:

    print('Warning: Record file overriden by -ro flag')

    ON_version = 'Overridden'
    ON_time1 = 'Overridden'
    ON_seq_file = 'Overridden'
    ON_seq_add_file = 'Overridden'
    ON_out = 'Overridden'
    ON_screen_type = 'Overridden'

    OFF_version = 'Overridden'
    OFF_time1 = 'Overridden'
    OFF_seq_file = 'Overridden'
    OFF_seq_add_file = 'Overridden'
    OFF_out = 'Overridden'
    OFF_screen_type = 'Overridden'


###############################################################################
# Pulls in ON and OFF counts and filters by defined threshold

print('Filtering reads')

# Retrieves filtered counts for auxilary function
ON, OFF, stats, time_zero = filterCounts(args.ON_file,
                                                args.OFF_file, args.thresh,
                                                args.zero_files, args.exclude)

# Outputs statistics
belowOFF, belowON, removed = stats

print('Total ON counts: ' + str(sum(ON.values())))
print('Total OFF counts: ' + str(sum(OFF.values())))
print('Missing/low ON elements: ' + str(belowON))
print('Missing/low OFF elements: ' + str(belowOFF))
print('Elements removed: ' + str(removed))
print('Number of distinct elements: ' + str(len(ON)))


###############################################################################
# Calculates enrichment values

print('Calculating enrichment values')

# Retrieves enrichment values from auxilary function
element_rhos, gene_rhos, neg_rhos, tar_rhos, gene_ref = enrich_all(ON,
		OFF, args.neg_name, args.split_mark, args.K, time_zero, args.back, args.bound_on, args.unbound_on)

print('Number of negative controls = ' + str(len(neg_rhos)))

if args.back == 'neg' and len(neg_rhos) == 0:
    sys.exit('No negative contols found.\n' + 
                'Change negative indicator with -n or --negative')

###############################################################################
# Writes output in human readable format

print('Outputing file')

with open(file_out + '_rhos.csv', 'w') as out_open:
    out_csv = csv.writer(out_open, delimiter=',', lineterminator='\n')

    for element, vals in element_rhos.items():

        # Writes to file
        out_csv.writerow([args.add + element, vals[0], vals[1], vals[2], vals[3], vals[4]])


###############################################################################
