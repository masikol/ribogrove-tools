#!/usr/bin/env python3

# TODO: add description

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
    '-i',
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
https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/.
It is better to filter this file with `filter_refseq_catalog.py` before running current script.
""",
    required=True
)

parser.add_argument(
    '-b',
    '--acc-blacklist',
    help="""A TSV file listing RefSeq Accession numbers (without version) to be discarded.
    The file should contain a header.
    Also, it should contains at least one column (of accession numbers). The second column (reason for rejection) is optional.
""",
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help="""a filtered TSV file (with header) with
    GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

args = parser.parse_args()


infpath = os.path.realpath(args.gi_2_acc_file)
filtered_catalog_fpath = os.path.realpath(args.refseq_catalog)
blacklist_fpath = os.path.realpath(args.acc_blacklist)
outfpath = os.path.realpath(args.outfile)


# Check existance of the input files

for fpath in (infpath, filtered_catalog_fpath, blacklist_fpath):
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

print(infpath)
print(filtered_catalog_fpath)
print(blacklist_fpath)
print()


# Read input
input_df = pd.read_csv(
    infpath,
    sep='\t'
)

# == Remove "WHOLE GENOME SHOTGUN SEQUENCE"s ==

print('Step 1')
print('Remove "whole genome shotgun" sequences')

def set_is_WGS(row):
    # Function sets `is_WGS` of a row to True if it's field `title` contains
    #    "WHOLE GENOME SHOTGUN SEQUENCE". Otherwise sets it to False
    row['is_WGS'] = 'WHOLE GENOME SHOTGUN SEQUENCE' in row['title'].upper()
    return row
# end def

print(f'number of rows before rm WGS = {input_df.shape[0]}')

input_df['is_WGS'] = np.repeat(None, input_df.shape[0])
input_df = input_df.apply(set_is_WGS, axis=1)
output_df = input_df[input_df['is_WGS'] == False]

print(f'number of rows after rm WGS = {output_df.shape[0]}')


# == Remove sequences added to RefSeq after the current release ==

print('\nStep 2')
print('Remove sequences added to RefSeq after the current release')

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



newly_added_df = output_df.query('not acc in @curr_release_accs')

newly_added_fpath = os.path.join(
    os.path.dirname(outfpath),
    'added_after_curr_release_' + os.path.basename(outfpath)
)

print(f'{newly_added_df.shape[0]} RefSeq records have been added to RefSeq since the current release')
print(f'Writing them to the file `{newly_added_fpath}`')

# Write records added after the current RefSeq release
newly_added_df.to_csv(
    newly_added_fpath,
    sep='\t',
    index=False,
    header=True,
    na_rep='NA',
    encoding='utf-8',
    columns=['gi_number', 'acc', 'title',]
)

# Filter remaining sequences
print('Filtering...')
output_df = output_df.query('acc in @curr_release_accs')
print(f'{output_df.shape[0]} RefSeq records remaining')


# == Remove sequences from the blacklist ==

print('\nStep 3')
print('Remove sequences from the blacklist')

backlist_df = pd.read_csv(
    blacklist_fpath,
    sep='\t'
)
blacklist_accs = set(backlist_df.iloc[:,0])

print('{} sequences found in the blacklist'.format(len(blacklist_accs)))
print('Filtering blacklist sequences (if there are any)')

def set_accs_no_version(row):
    row['acc_no_version'] = row['acc'].partition('.')[0]
    return row
# end def

output_df['acc_no_version'] = np.repeat(None, output_df.shape[0])
output_df = output_df.apply(set_accs_no_version, axis=1)
output_df = output_df.query('not acc_no_version in @blacklist_accs')

# Clear
output_df = output_df.drop(columns=['acc_no_version'])

print(f'{output_df.shape[0]} RefSeq records remaining')


# == Write output ==

output_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    na_rep='NA',
    encoding='utf-8',
    columns=['gi_number', 'acc', 'title',]
)

print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
