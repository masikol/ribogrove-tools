#!/usr/bin/env python3

# The script annotates gene sequences -- adds taxonomy and category info.
#
# Input files:
# 1. Input fasta file (-f/--fasta-seqs-file).
# 2. TSV file mapping genes seqIDs to taxonomy (-t/--per-gene-taxonomy-file)
# 3. TSV file mapping genes seqIDs to categories (-c/--categories-file)
#
# Output files:
# 1. Annotated fasta file (-o/--outfile)

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import argparse

import pandas as pd
from Bio import SeqIO


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--fasta-seqs-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-t',
    '--per-gene-taxonomy-file',
    help='TSV file (with header) containing per-gene taxonomy',
    required=True
)

parser.add_argument(
    '-c',
    '--categories-file',
    help='TSV file (with header) containing categories info',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='output file of annotated gene sequences',
    required=True
)

args = parser.parse_args()


# For convenience
in_fasta_fpath = os.path.abspath(args.fasta_seqs_file)
tax_fpath = os.path.abspath(args.per_gene_taxonomy_file)
cat_fpath = os.path.abspath(args.categories_file)
outfpath = os.path.abspath(args.outfile)


# Check existance of all input files and dependencies
for fpath in (in_fasta_fpath, tax_fpath, cat_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directory if needed
if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if


print(in_fasta_fpath)
print(tax_fpath)
print(cat_fpath)
print()


tax_sep = ';'

# == Proceed ==

# Read taxonomy file
tax_df = pd.read_csv(
    tax_fpath,
    sep='\t',
    dtype={
        'seqID': str,
        'ass_id': pd.Int32Dtype(),
        'accs': str,
        'taxID': pd.Int32Dtype(),
        'tax_name': str,
        'genus': str,
        'family': str,
        'order': str,
        'class': str,
        'phylum': str,
        'superkingdom': str,
    }
)

# Read categories file
cat_df = pd.read_csv(
    cat_fpath,
    sep='\t',
    dtype={
        'ass_id': pd.Int32Dtype(),
        'seqID': str,
        'category': pd.Int8Dtype(),
    }
)

# Count sequences
n_seqs = len(tuple(SeqIO.parse(in_fasta_fpath, 'fasta')))

# In order to select sequences quickly
tax_df.index = tax_df['seqID']


with open(outfpath, 'wt') as outfile:

    seq_records = SeqIO.parse(in_fasta_fpath, 'fasta')

    for i, seq_record in enumerate(seq_records):

        print(f'\r Doing {i+1}/{n_seqs}: {seq_record.id}', end=' '*10)

        # Select line of taxonomy DF for current sequence
        curr_tax_record = tax_df.loc[seq_record.id, ]

        # Form taxonomy (lineage) string
        taxonomy = tax_sep.join(
            (
                'NA' if pd.isnull(curr_tax_record['superkingdom']) else curr_tax_record['superkingdom'],
                'NA' if pd.isnull(curr_tax_record['phylum']) else curr_tax_record['phylum'],
                'NA' if pd.isnull(curr_tax_record['class']) else curr_tax_record['class'],
                'NA' if pd.isnull(curr_tax_record['order']) else curr_tax_record['order'],
                'NA' if pd.isnull(curr_tax_record['family']) else curr_tax_record['family'],
                'NA' if pd.isnull(curr_tax_record['genus']) else curr_tax_record['genus'],
            )
        )

        # Form taxonomy name
        if not pd.isnull(curr_tax_record['tax_name']):
            tax_name = curr_tax_record['tax_name'].replace(' ', '_')
        else:
            tax_name = 'no_taxonomy_name'
        # end if

        # Select category
        category = cat_df[cat_df['seqID'] == seq_record.id]['category'].values[0]

        # Form header for output fasta file
        seq_record.description = f'{seq_record.id} {tax_name} {tax_sep}{taxonomy}{tax_sep} category:{category}'

        outfile.write(f'>{seq_record.description}\n{seq_record.seq}\n')
    # end for
# end with

print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
