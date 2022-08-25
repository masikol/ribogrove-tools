#!/usr/bin/env python3

# The script merges two files together:
# 1) a TSV file, which is output of the script `assembly2gi_numbers.py` (`-c` option);
# 2) a TSV file, which is output of the script `gis_to_accs.py` (`-c` option).

# The output file (`-o` option) is a TSV file of 4 columns (`ass_id`, `gi_number`, `acc`, `title`).
#   In this file, every line corresponds to a single RefSeq record, and the line
#   contains Assembly ID (`ass_id`), RefSeq GI number (`gi_number`), ACCESSION.VERSION (`acc`),
#   and the title of the record (`title`).

## Command line arguments
### Input files:
# 1. `-s / --assm-2-gi-file` -- input TSV file mapping Assembly IDs to RefSeq GI numbers. Mandatory.
# 2. `-c / --gi-2-acc-file` -- input TSV file mapping RefSeq GI numbers to RefSeq ACCESSION.VERSIONs
#   and titles. Mandatory.
# 3. `-a / --refseq-catalog` -- A RefSeq "catalog" file of the current release.
#   This is the file `RefSeq-releaseXXX.catalog.gz` from here:
#   https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/.
#   It is better to filter this file with `filter_refseq_catalog.py` before running current script.

### Output files:
# 1. `-o / --outfile` -- output TSV files where Assembly IDs are mapped to ACCESSION.VERSIONs
#   and RefSeq Titles. Mandatory.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import gzip
import argparse

import numpy as np
import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-s',
    '--assm-2-gi-file',
    help='TSV file (with header) with Assembly IDs and GI numbers separated by tabs',
    required=True
)

parser.add_argument(
    '-c',
    '--gi-2-acc-file',
    help="""TSV file (with header) with
    GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='file mapping Assembly IDs to RefSeq accession and titles',
    required=True
)

args = parser.parse_args()


assm_2_gi_fpath = os.path.realpath(args.assm_2_gi_file)
gi_2_acc_fpath = os.path.realpath(args.gi_2_acc_file)
outfpath = os.path.realpath(args.outfile)


# Check existance of the input files

for fpath in (assm_2_gi_fpath, gi_2_acc_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# end for

if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

print(assm_2_gi_fpath)
print(gi_2_acc_fpath)
print()


# Read input
ass_2_gi_df = pd.read_csv(
    assm_2_gi_fpath,
    sep='\t'
)
print('Assembly IDs and GI numbers:')
print('{} rows'.format(ass_2_gi_df.shape[0]))
print(ass_2_gi_df.head())
print()

gi_2_acc_df = pd.read_csv(
    gi_2_acc_fpath,
    sep='\t'
)
print("GI numbers, 'ACCESION.VERSION's, and titles:")
print('{} rows'.format(gi_2_acc_df.shape[0]))
print(gi_2_acc_df.head())
print()


# == Merge Assembly IDs with RefSeq ACCESSION.VERSION ==

merged_df = ass_2_gi_df.merge(gi_2_acc_df, on='gi_number', how='right')

print('Merged dataframe:')
print("Assembly IDs, GI numbers, 'ACCESION.VERSION's, and titles:")
print('{} rows'.format(merged_df.shape[0]))
print(merged_df.head())


# == Write output ==

merged_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    na_rep='NA',
    encoding='utf-8',
    columns=['ass_id', 'gi_number', 'acc', 'title',]
)

print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
