###############################################################################
# Josh Tycko
# 03/24/2019
# Adapted from casTLE by David Morgens
###############################################################################
# Import neccessary modules

'''
Function for aligning fastq files.
'''

import subprocess
import shlex
import sys
import csv
import os
import argparse
import HTSeq
from collections import defaultdict
import time


############################################################################### 
# Version number

current_version = '0.5'


###############################################################################
# Retrieve alignment indexes by screen type

name2index = {}

# Reads in index file names
index_file = os.path.join('Indices', 'screen_type_index.txt')

with open(index_file, 'r') as index_open:
    index_csv = csv.reader(index_open, delimiter='\t')
    for row in index_csv:
        name2index[row[0]] = row[1]


###############################################################################
# Parses input using argparse module

# Initiates argument parser
parser = argparse.ArgumentParser(description='Align FASTQ and make count file')

# Non-optional arguments:
parser.add_argument('file_in', help='File base for input fastq files',
                    type=str)

parser.add_argument('name', help='Name for output files',
                    type=str)

parser.add_argument('screen_type', help='The screen type', type=str,
                        choices=sorted(name2index.keys()))

# Optional arguments:
parser.add_argument('-m', '--mismatch', dest='mismatch',
                    help='The number of tolerated mismatches; default is 0',
                    type=str, default='0')

parser.add_argument('-t', '--trim', dest='trim',
                    help='Select the number of bases to trim off start of read; default is 25',
                    type=int, default=25)

parser.add_argument('-l', '--length', dest='read_length',
                    help='Select the number of bases to align; default is 150',
                    type=int, default=150)

parser.add_argument('-fi', '--filter',
                    help='Flag to filter too short reads', action='store_true')

parser.add_argument('-b', '--bowtie',
                    help='Location of Bowtie aligner; default is bowtie', type=str,
                    default='bowtie')

parser.add_argument('-a', '--add_file',
          help='Location of additional FASTQ files, if any')

parser.add_argument('-s', '--strand', type=str, choices=['-','+'],
          help='Filters reads by alignment strand; default is +', default='+')

parser.add_argument('-p', '--process',
                help='Number of processors to use; default is 12', type=str,
                default="12")

# Saves input to args object
args = parser.parse_args()


###############################################################################
# Makes checks for failed arguments

# Determines output file
file_out = os.path.join('Data', args.name)

try:
    with open(file_out + '_record.txt','w') as file_open:
        pass
except:
    sys.exit('Cannot write to output file: \n' + file_out + '\n'
                + 'Use -of or --override_file to change')

# Determines if bowtie location is valid
if not os.path.exists(args.bowtie):
    found = False

    # Looks for shortcuts i.e. bowtie
    for path in os.environ["PATH"].split(os.pathsep):
        path = path.strip('"')
        exe = os.path.join(path, args.bowtie)
        if os.path.exists(exe):
            found = True
            break

    if not found:
        sys.exit('Bowtie aligner not found at: \n' + args.bowtie + '\n'
                    + 'Use -b or --bowtie to change')

# Locates index files
index = name2index[args.screen_type]


###############################################################################
# Gathers reads using unix shell commands

# Calls the unix shell to gather all read files starting with file base
print('Finding reads')

# If additional files, collect these too
if args.add_file:
    try:
        subprocess.check_call('cat ' + args.file_in + '* ' + args.add_file
                                + '*' + ' > ' + file_out + '_all.fastq.gz',
                                shell=True)
    except:
        sys.exit('Shell error')
else:
    try:
        subprocess.check_call('cat ' + args.file_in + '*'+' > '
                                + file_out + '_all.fastq.gz', shell=True)
    except:
        sys.exit('Shell error')


# Calls the unix shell to unzip newly aggregated read file
print('Unzipping reads')

try:
    subprocess.check_call('gunzip ' + file_out + '_all.fastq.gz -f', shell=True)
except:
    sys.exit('Shell error')


###############################################################################
# Trims reads using HTSeq

total_reads = 0  # Stores count of number of reads in file
warn = True
short = 0

# Moves over all reads, trims off the primer sequence at start of read, then further trims them down to the read length (args.read_length).
# Then writes trimmed reads to new fastq file

print('Trimming reads')

with open(file_out + '_trimmed.fastq', 'w') as align_file:

    for read in HTSeq.FastqReader(file_out + '_all.fastq'):
        total_reads += 1
        trimmed_read = read[args.trim : args.trim + args.read_length]

        # Gives a warning if read shorter than trimmed length
        if len(trimmed_read) != args.read_length:
            short += 1

            if warn:
                print('Warning: Read shorter than indicated length\n'
                            + 'Read length: ' + str(len(trimmed_read)))
                warn = False

            if args.filter:
                continue

        trimmed_read.write_to_fastq_file(align_file)

if not warn:
    print(str(short) + ' reads too short.\n'
                + 'Use -f or --filter to remove')


###############################################################################
# Runs the alignment with bowtie

# Create files for input, output, and index for bowtie
fastqFile = file_out + '_trimmed.fastq'
mapFile = file_out + '.map'
unmapFile = file_out + '.unmapped'

print('Mapping reads')

# Calls bowtie in shell.
# -v is the number of tolerate mismatches.
# -p is the number of processors to use.
# -q is the input file and indicates that it is fastq.
# mapFile stores the output alignments.
# --un stores the reads with no alignments
subprocess.call(shlex.split(" ".join([args.bowtie,"-v", args.mismatch,
                            "-a -p", args.process, index, "-q", fastqFile,
                            mapFile, "--un", unmapFile])))


###############################################################################
# Parses bowtie output

# Stores counts of reads by entry
counts = defaultdict(int)

print('Analyzing reads')

# Opens bowtie output, checks to make sure correct strand 
# And then counts reads by element.

total_counts = 0
ambig_counts = 0
last_read = ''
ambig = 0
warn = True

with open(mapFile, 'r') as map_open:
    map_csv = csv.reader(map_open, delimiter='\t')

    for line in map_csv:

        name = line[2]
        strand = line[1]
        read = line[0]

        # Checks if alignment is to the correct strand
        if args.strand and strand != args.strand:
            continue

        counts[name] += 1
        total_counts += 1

        # Counts multimapped reads
        if read == last_read:
            ambig_counts += 1
            ambig += 1
        else:
            ambig = 0

        last_read = read

        # Warning for extreme multimapping reads
        if ambig > 100 and warn:
            print('Warning: Possible alignment to constant region')
            warn = False

ambig_pct = round(float(ambig_counts)/total_counts*100.0,3)
print(str(total_counts) + ' total counts')
print(str(ambig_counts) + ' ambiguous counts')
print(str(ambig_pct) + ' ambiguous percent')


###############################################################################
# Write output and records

# Writes output as tab deliminated file
with open(file_out + '_counts.csv', 'w') as out_open:
    out_csv = csv.writer(out_open, delimiter=',')
    for name in sorted(counts):
        out_csv.writerow([name, counts[name]])

# Stores record file for downstream analysis
with open(file_out + '_record.txt', 'w') as rec_open:
    rec_csv = csv.writer(rec_open, delimiter='\t')
    rec_csv.writerow(['makeCounts.py', current_version])
    rec_csv.writerow(['Date', time.strftime("%d:%m:%Y")])
    rec_csv.writerow(['Sequencing files', args.file_in])
    rec_csv.writerow(['Additional Files', args.add_file])
    rec_csv.writerow(['Output File', file_out])
    rec_csv.writerow(['Screen Type', args.screen_type])
    rec_csv.writerow(['Mismatches allowed', args.mismatch])
    rec_csv.writerow(['Length of alignment', args.read_length])
    rec_csv.writerow(['Trim from start', args.trim])
    rec_csv.writerow(['Reads input', total_reads])
    rec_csv.writerow(['Total aligned reads', total_counts])
    rec_csv.writerow(['Ambiguous aligned reads', ambig_counts])
    rec_csv.writerow(['Alignment %', round(float(total_counts-ambig_counts)/total_reads*100.0,3)])
    rec_csv.writerow(['Ambiguous %', ambig_pct])
    

# Stores record file for permanent record
with open(os.path.join('Records', 'makeCounts' + time.strftime("%Y%m%d%H%M%S")), 'w') as back_open:
    rec_csv = csv.writer(back_open, delimiter='\t')
    rec_csv.writerow(['makeCounts.py', current_version])
    rec_csv.writerow(['Date', time.strftime("%d:%m:%Y")])
    rec_csv.writerow(['Sequencing files', args.file_in])
    rec_csv.writerow(['Additional Files', args.add_file])
    rec_csv.writerow(['Output File', file_out])
    rec_csv.writerow(['Screen Type', args.screen_type])
    rec_csv.writerow(['Mismatches allowed', args.mismatch])
    rec_csv.writerow(['Length of alignment', args.read_length])
    rec_csv.writerow(['Trim from start', args.trim])
    rec_csv.writerow(['Reads input', total_reads])
    rec_csv.writerow(['Total aligned reads', total_counts])
    rec_csv.writerow(['Ambiguous aligned reads', ambig_counts])
    rec_csv.writerow(['Alignment %', round(float(total_counts-ambig_counts)/total_reads*100.0,3)])
    rec_csv.writerow(['Ambiguous %', ambig_pct])


###############################################################################
# Compresses and/or erases files using linux shell commands

print('Deleting gathered reads')
try:
    subprocess.check_call('rm ' + file_out + '_all.fastq', shell=True)
except:
    sys.exit('Shell error')
try:
    subprocess.check_call('rm ' + file_out + '_trimmed.fastq', shell=True)
except:
    sys.exit('Shell error')

print('Compressing mapped reads')
try:
    subprocess.check_call('gzip ' + file_out + '.map -f', shell=True)
except:
    sys.exit('Shell error')

print('Compressing unmapped reads')
try:
    subprocess.check_call('gzip ' + file_out + '.unmapped -f', shell=True)
except:
    sys.exit('Shell error')

print('Finished')


###############################################################################
