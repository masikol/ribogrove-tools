#!/usr/bin/env python3

# The script merges all info about the gene sequences (nucleotide composition, category, taxonomy)
#   into a single TSV file.

## Command line arguments
### Input files:
# 1. `-b / --bases-file` -- a input TSV file with counted bases ouputed
#   by the script `count_bases.py`.
#   This is the file, which is the output of th script `count_bases.py`. Mandatory.
# 2. `-c / --categories-file` -- the per-gene TSV file of categories info.
#   This is the file, which is the output of th script `assign_genome_categories.py`. Mandatory.
# 3. `-t / --taxonomy-file` -- the per-gene taxonomy file.
#   This is the file, which is the output of th script `add_taxonomy_names.py`. Mandatory.

### Output files:
# 1. `-o / --outfile` -- an output TSV file. Mandatory.


import os
import time

def get_time():
    return time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())
# end def

print(
    '\n|=== {} STARTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)


# == Parse arguments ==
import argparse

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


# == Import them now ==
import sys

import numpy as np
import pandas as pd

from src.ribogrove_seqID import parse_asm_acc


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


def set_asm_acc(row):
    row['asm_acc'] = parse_asm_acc(row['seqID'])
    return row
# end def


# seqID   a       t       g       c  [...]     len
bases_df = pd.read_csv(bases_fpath, sep='\t')
bases_df['asm_acc'] = np.repeat('', bases_df.shape[0])
bases_df = bases_df.apply(set_asm_acc, axis=1)

# asm_acc category
categories_df = pd.read_csv(categories_fpath, sep='\t')

# seqID   asm_acc  accs    taxID   tax_name        genus   family  order   class   phylum  domain
taxonomy_df = pd.read_csv(taxonomy_fpath, sep='\t')


print('Merging...')
merged_df = bases_df \
    .merge(categories_df[['asm_acc', 'category']], on='asm_acc', how='left') \
    .merge(taxonomy_df[['asm_acc', 'taxid', 'organism_name', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']],
        on='asm_acc',
        how='left'
    )
print('  done!\n')

# Reorder columns
merged_df = merged_df[
    [
        'seqID',
        'asm_acc', 'category',
        'taxid', 'organism_name', 'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain',
        'len', 'a', 't', 'g', 'c', 'r', 'y', 'w', 's', 'k', 'm', 'h', 'v', 'b', 'd', 'n',
    ]
]

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
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
