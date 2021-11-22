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

parser.add_argument(
    '-a',
    '--refseq-catalog',
    help="""A RefSeq "catalog" file of the current release.
This is the file `RefSeq-releaseXXX.catalog.gz` from here:
https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/""",
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
filtered_catalog_fpath = os.path.realpath(args.refseq_catalog)
outfpath = os.path.realpath(args.outfile)


# Check existance of the input files

for fpath in (assm_2_gi_fpath, gi_2_acc_fpath, filtered_catalog_fpath):
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
print(filtered_catalog_fpath)
print()


# Read input
ass_2_gi_df = pd.read_csv(
    assm_2_gi_fpath,
    sep='\t'
)

gi_2_acc_df = pd.read_csv(
    gi_2_acc_fpath,
    sep='\t'
)


print(ass_2_gi_df.shape)
print(ass_2_gi_df.head())

print(gi_2_acc_df.shape)
print(gi_2_acc_df.head())


def set_is_WGS(row):
    # Function sets `is_WGS` of a row to True if it's field `title` contains
    #    "WHOLE GENOME SHOTGUN SEQUENCE". Otherwise sets it to False
    row['is_WGS'] = 'WHOLE GENOME SHOTGUN SEQUENCE' in row['title'].upper()
    return row
# end def set_is_WGS


# == Remove "WHOLE GENOME SHOTGUN SEQUENCE"s ==

print(f'number of rows before rm WGS = {gi_2_acc_df.shape[0]}')

gi_2_acc_df['is_WGS'] = np.repeat(None, gi_2_acc_df.shape[0])
gi_2_acc_df = gi_2_acc_df.apply(set_is_WGS, axis=1)
gi_2_acc_df = gi_2_acc_df[gi_2_acc_df['is_WGS'] == False]

print(f'number of rows after rm WGS = {gi_2_acc_df.shape[0]}')


# == Merge Assembly IDs with RefSeq ACCESSION.VERSION ==

merged_df = ass_2_gi_df.merge(gi_2_acc_df, on='gi_number', how='right')

print('MERGED DATAFRAME:')
print(merged_df.shape)
print(merged_df.head())


print('\nRemoving sequences added to RefSeq after the current release')

# Read the catalog file
if filtered_catalog_fpath.endswith('.gz'):
    open_func = gzip.open
else:
    open_func = open
# end if

curr_release_accs = set()

with open_func(filtered_catalog_fpath, 'rt') as catalog_file:
    acc_column_index = 2
    dir_column_index = 3
    separator = '\t'

    for line in catalog_file:
        line_vals = line.split(separator)
        curr_release_accs.add(line_vals[acc_column_index])
    # end for
# end with
print('done\n')


newly_added_df = merged_df.query('not acc in @curr_release_accs')

newly_added_fpath = os.path.join(
    os.path.dirname(outfpath),
    'added_after_curr_release_' + os.path.basename(outfpath)
)

print(f'{newly_added_df.shape[0]} RefSeq records have been added to RefSeq since the current release')
print(f'Writing them to the file `{newly_added_fpath}`')


merged_df = merged_df.query('acc in @curr_release_accs')
print(f'{merged_df.shape[0]} RefSeq records remaining')
print(f'Writing them to the output file `{outfpath}`')

# Write output
merged_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    na_rep='NA',
    encoding='utf-8',
    columns=['ass_id', 'gi_number', 'acc', 'title',]
)


# Write records added after the current RefSeq release
newly_added_df.to_csv(
    newly_added_fpath,
    sep='\t',
    index=False,
    header=True,
    na_rep='NA',
    encoding='utf-8',
    columns=['ass_id', 'gi_number', 'acc', 'title',]
)

print('\nCompleted!')
print(outfpath)
print(newly_added_fpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
