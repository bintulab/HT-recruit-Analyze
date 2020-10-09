###############################################################################
# Josh Tycko
# 03/25/2019
# Adapted from casTLE (David Morgens)
################################################################################
# Import neccessary modules

import argparse
import sys
import csv
import subprocess
import shlex
import os


##############################################################################
# Initiates argument parser

parser = argparse.ArgumentParser('Make index for alignment')

parser.add_argument('oligo_file', help='Input oligo file; csv format', type=str)

parser.add_argument('short_name', help='The screen type for reference', type=str)

parser.add_argument('full_name', help='Name for output files', type=str)

# Optional arguments: base trimming of fasta

parser.add_argument('-s','--strim', help='Trim bases from start; default is 30',
                        default=30, type=int)

parser.add_argument('-e', '--etrim', help='Trim bases from end; default is 30',
                        default=30, type=int)

parser.add_argument('-n', '--nums', help='Number of oligos to output',
                        default=3, type=int)

parser.add_argument('-o', '--override', help='Flag to override existing indexes',
                        action='store_true')

parser.add_argument('-t', '--test', help="Flag to not run bowtie",
                        action='store_false')

parser.add_argument('-b', '--bowtie',
                    help='Location of Bowtie aligner; default is /usr/bin/bowtie-build', type=str,
                    default='/usr/bin/bowtie-build')

args = parser.parse_args()


##############################################################################
#

index_file = os.path.join('Indices', 'screen_type_index.txt')
index = []

# Check whether screen type and name already exist
with open(index_file, 'rU') as index_open:

    index_csv = csv.reader(index_open, delimiter='\t')

    for line in index_csv:

        name, location = line

        if name == args.short_name or location == args.full_name:

            if not args.override:
                sys.exit('Error: Screen name taken')

            else:
                print('Warning: Previous screen overwritten')

        else:
            index.append((name, location))


##############################################################################
# Convert oligo file into bowtie-compatible fasta

with open(args.oligo_file, 'r') as oligo_file:

    oligo_csv = csv.reader(oligo_file, delimiter=',')
    oligo_list = []
    
    if args.etrim > 0:
        for line in oligo_csv:
            oligo = line[1][args.strim: -args.etrim]
            oligo_list.append(['>' + line[0], oligo])

    else:
        for line in oligo_csv:
            oligo = line[1][args.strim: ]
            oligo_list.append(['>' + line[0], oligo])

if args.test:

    for t in range(args.nums):
        oligo = oligo_list[t][1]
        print('Sample oligo: ' + oligo)

    for t in range(args.nums):
        oligo = oligo_list[t][1]
        print('Sample length: ' + str(len(oligo)))

    sys.exit('Warning: Files not created.\n' +
                'Use -t or --test to create files')


##############################################################################
# 

fasta_location = os.path.join('Indices', 'temp_fasta.fna')

with open(fasta_location, 'w') as fasta_open:

    fasta_csv = csv.writer(fasta_open, delimiter='\n')

    for i in oligo_list:
        fasta_csv.writerow(i)


##############################################################################
# Call bowtie-build to build new index

try:
    subprocess.check_call(args.bowtie + ' ' + fasta_location + ' ' + 
                                            os.path.join('Indices', args.full_name),
                                                                    shell=True)
except:
    sys.exit()

# delete fasta files
try:
    subprocess.check_call('rm ' + fasta_location, shell=True)
except:
    sys.exit()


##############################################################################
# Add new index to screentype_index file

index_name = os.path.join('Indices', args.full_name)
index.append([args.short_name, index_name])
index.sort(key=lambda x: x[0])

with open(index_file, 'w') as index_open:

    index_csv = csv.writer(index_open, delimiter='\t')

    for ind in index:
        index_csv.writerow(ind)


##############################################################################
