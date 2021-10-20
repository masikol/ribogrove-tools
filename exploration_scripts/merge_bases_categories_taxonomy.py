#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# The script merges all info about gene sequences (nucleotide composition, category, taxonomy)
#   into a single TSV file

# Input files
# 1. -b/--bases-file -- input TSV file with counted bases ouputed by the script count_bases.py
# 2. -c/--categories-file -- per-gene TSV file of categories info
# 3. -t/--taxonomy-file -- per-gene taxonomy file

# Output files
# -o/--outfile -- output TSV file

import os
import sys
import argparse

import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-b',
    '--bases-file',
    help='input TSV file with counted bases ouputed by the script count_bases.py',
    required=True
)

parser.add_argument(
    '-c',
    '--categories-file',
    help='per-gene TSV file of categories info',
    required=True
)

parser.add_argument(
    '-t',
    '--taxonomy-file',
    help='per-gene taxonomy file',
    required=True
)


# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='output merged TSV file',
    required=True
)

args = parser.parse_args()

# For convenience
bases_fpath = os.path.abspath(args.bases_file)
categories_fpath = os.path.abspath(args.categories_file)
taxonomy_fpath = os.path.abspath(args.taxonomy_file)
outfpath = os.path.abspath(args.outfile)

for f in [bases_fpath, categories_fpath, taxonomy_fpath]:
    if not os.path.exists(f):
        print(f'Error: input file `{f}` does not exist.')
        sys.exit(1)
    # end if
# end for

if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Cannot create output directory `{os.path.dirname(outfpath)}`')
        print(str(err))
        sys.exit(1)
    # end try
# end if

print(bases_fpath)
print(categories_fpath)
print(taxonomy_fpath)
print()


# seqID   a       t       g       c       len
bases_df = pd.read_csv(bases_fpath, sep='\t')

# ass_id  seqID   category
categories_df = pd.read_csv(categories_fpath, sep='\t')

# seqID   ass_id  accs    taxID   tax_name        genus   family  order   class   phylum  superkingdom
taxonomy_df = pd.read_csv(taxonomy_fpath, sep='\t')


print('Merging...')
merged_df = bases_df \
    .merge(categories_df[['seqID', 'category', 'ass_id']], on='seqID', how='left') \
    .merge(taxonomy_df[['seqID', 'taxID', 'tax_name', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'superkingdom']],
        on='seqID',
        how='left'
    )
print('  done!\n')

print('Merged dataframe:')
print(merged_df.shape)
print(merged_df.head())

merged_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    encoding='utf-8',
    na_rep='NA'
)

print('Completed!')
print(outfpath)
