###############################################################################
# Josh Tycko
# 03/25/2019
###############################################################################
# Import neccessary modules

'''
Combine rhos and make figures 
'''

from __future__ import division
import csv
import time
import argparse
from screenFun import *
import sys
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['font.sans-serif'] = "Arial"
mpl.rcParams['font.family'] = "sans-serif"
mpl.rcParams['pdf.fonttype'] = 42

###############################################################################    
# Version number

current_version = '1.0'


###############################################################################    

# Parses input using argparse module

# Initiates input parser
parser = argparse.ArgumentParser(description='Combines rhos (enrichment ratios) from replicates')

# Non-optional arguments:
parser.add_argument('R1_file', help='File for R1 rhos', type=str)

parser.add_argument('R2_file', help='File for R2 rhos', type=str)

parser.add_argument('name', help='Name for output files', type=str)

# Optional arguments
parser.add_argument('-of', '--override_file', action='store_true',
                help='Override Result file output location.')

# parser.add_argument('-n', '--neg_name',
#                 help='Symbol used to denote negative controls',
#                 type=str, default='0')

# parser.add_argument('-s', '--split', dest='split_mark', help='Delimiter for element name',
#                 type=str, default='_')

# parser.add_argument('-x', '--exclude', type=str, help='Excludes substrings', nargs='+')


# Saves all input to object args
args = parser.parse_args()


###############################################################################
# Processes and checks input arguments

# Creates output file name
if args.override_file:
    file_out = args.name
    figure_out = args.name
else:
    file_out = os.path.join('Results', args.name)
    figure_out = os.path.join('Results/Figures', args.name)

# Checks if can write to output
try:
    with open(file_out + '_record.txt', 'w') as out_open:
        pass

    os.remove(file_out + '_record.txt')

except:
    sys.exit('Cannot write to output file:\n' + file_out + '\n'
                + 'Use -of or --override_file to change')


###############################################################################
# Load and combine rhos from individual replicates

df1 = pd.read_csv(args.R1_file, names = ['label', 'R1', 'countsOFF_R1', 'countsON_R1', 'Poff R1', 'Poff_corrected R1'], header = None)
df2 = pd.read_csv(args.R2_file, names = ['label', 'R2', 'countsOFF_R2', 'countsON_R2', 'Poff R2', 'Poff_corrected R2'], header = None)

df = pd.merge(df1, df2, how = 'outer', on = 'label')

df['Avg'] = df[['R1','R2']].mean(axis = 1)
df['Standard Error'] = df[['R1','R2']].sem(axis = 1)

###############################################################################
# Writes output in human readable format

df.to_csv(file_out + '_combo.csv', index=None, header=True)

###############################################################################
# Make reproducibility figures
g = sns.jointplot(x ='R1', y = 'R2', data = df, alpha = 0.3, linewidth = 1, space = 0)
g.set_axis_labels('Replicate 1 (log2 OFF:ON)', 'Replicate 2 (log2 OFF:ON)')
plt.savefig(figure_out+ '_reps_rhos.pdf', transparent = True)
plt.close()

g = sns.jointplot(x ='countsOFF_R1', y = 'countsOFF_R2', data = df, alpha = 0.3, linewidth = 1, space = 0)
g.set_axis_labels('Replicate 1 (counts OFF)', 'Replicate 2 (counts OFF)')
plt.savefig(figure_out+ '_reps_OFFcounts.pdf', transparent = True)
plt.close()

g = sns.jointplot(x ='countsON_R1', y = 'countsON_R2', data = df, alpha = 0.3, linewidth = 1, space = 0)
g.set_axis_labels('Replicate 1 (counts ON)', 'Replicate 2 (counts ON)')
plt.savefig(figure_out+ '_reps_ONcounts.pdf', transparent = True)
plt.close()
