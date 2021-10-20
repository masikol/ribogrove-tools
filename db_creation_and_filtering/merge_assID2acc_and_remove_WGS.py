#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# The script merges TSV file, which is output of the script assembly2refseq_id.py (-s/--assm-2-gi-file)
#   and TSV file, which is output of the script gis_to_accs.py (-c/--gi-2-acc-file) on column `refseq_id`.
# Output (file -o/--outfile) is a TSV file of 4 columns (ass_id, refseq_id, acc, title).

# Input files
# 1. -s/--assm-2-gi-file -- input TSV file mapping Assembly IDs to RefSeq GI numbers.
# 2. -c/--gi-2-acc-file -- input TSV file mapping RefSeq GI numbers to RefSeq ACCESSION.VERSIONs.

# Output files
# 1. -o/--outfile -- output TSV files where Assembly IDs are mapped to ACCESSION.VERSIONs and RefSeq Titles.


import os
import sys
import argparse

import numpy as np
import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

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
    '-o',
    '--outfile',
    help='file mapping Assembly IDs to RefSeq accession and titles',
    required=True
)

args = parser.parse_args()


assm_2_gi_fpath = os.path.realpath(args.assm_2_gi_file)
gi_2_acc_fpath = os.path.realpath(args.gi_2_acc_file)
outfpath = os.path.realpath(args.outfile)


# Check existance of input file -s/--assm-2-gi-file
if not os.path.exists(assm_2_gi_fpath):
    print(f'Error: file `{assm_2_gi_fpath}` does not exist!')
    sys.exit(1)
# end if

# Check existance of input file -c/--gi-2-acc-file
if not os.path.exists(gi_2_acc_fpath):
    print(f'Error: file `{gi_2_acc_fpath}` does not exist!')
    sys.exit(1)
# end if

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

merged_df = ass_2_gi_df.merge(gi_2_acc_df, on='refseq_id', how='right')

print('MERGED: DATAFRAME')
print(merged_df.shape)
print(merged_df.head())

# Write output
merged_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    na_rep='NA',
    encoding='utf-8',
    columns=['ass_id', 'refseq_id', 'acc', 'title',]
)

print('\nCompleted!')
print(outfpath)
