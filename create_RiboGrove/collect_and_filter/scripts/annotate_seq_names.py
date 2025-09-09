#!/usr/bin/env python3

# The script annotates gene sequences -- adds taxonomy and category info.

## Command line arguments

### Input files:
# 1. `-f / --fasta-seqs-file` -- input fasta file to be annotated.
#   This file is the output of the script `make_final_seqs.py`. Mandatory.
# 2. `-t / --taxonomy` -- a TSV file mapping RiboGrove seqIDs to taxonomy.
#   This file is the output of the script `add_taxonomy_names.py`. Mandatory.
# 3. `-c / --categories-file` -- a TSV file mapping RiboGrove seqIDs to categories.
#   This file is the output of the script `assign_genome_categories.py`. Mandatory.

### Output files:
# 1. `-o / --outfile` -- an annotated ouput fasta file. Mandatory.


import os
from src.rg_tools_time import get_time

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
    '-f',
    '--fasta-seqs-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-t',
    '--taxonomy-file',
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


# == Import them now ==
import sys

import pandas as pd
from Bio import SeqIO

from src.ribogrove_seqID import parse_asm_acc


# For convenience
in_fasta_fpath = os.path.abspath(args.fasta_seqs_file)
tax_fpath = os.path.abspath(args.taxonomy_file)
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


TAX_SEP = ';'


def make_species_name(raw_species):
    species_name = raw_species
    
    if raw_species.startswith('Candidatus '):
        species_name = species_name.replace('Candidatus ', '')
    # end if

    species_name = species_name.partition(' ')[2]

    return species_name
# end def


def get_taxonomy(asm_acc, tax_df):
    global TAX_SEP

    # Select line of taxonomy DF for current sequence
    curr_tax_record = tax_df.loc[asm_acc,]

    raw_species_name = 'NA' if pd.isnull(curr_tax_record['Species']) else curr_tax_record['Species']
    species_name = make_species_name(raw_species_name)

    tax_names = [
        'NA' if pd.isnull(curr_tax_record[ 'Domain']) else curr_tax_record[ 'Domain'],
        'NA' if pd.isnull(curr_tax_record['Kingdom']) else curr_tax_record['Kingdom'],
        'NA' if pd.isnull(curr_tax_record[ 'Phylum']) else curr_tax_record[ 'Phylum'],
        'NA' if pd.isnull(curr_tax_record[  'Class']) else curr_tax_record[  'Class'],
        'NA' if pd.isnull(curr_tax_record[  'Order']) else curr_tax_record[  'Order'],
        'NA' if pd.isnull(curr_tax_record[ 'Family']) else curr_tax_record[ 'Family'],
        'NA' if pd.isnull(curr_tax_record[  'Genus']) else curr_tax_record[  'Genus'],
        species_name,
    ]

    prefixes = ['d__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__',]

    tax_names_with_prefixes = [
        '{}{}'.format(prefix, name) for prefix, name in zip(prefixes, tax_names)
    ]

    # Form taxonomy (lineage) string
    taxonomy_str = TAX_SEP.join(tax_names_with_prefixes)
    taxonomy_str = taxonomy_str.replace(' ', '_')

    return taxonomy_str
# end def


def get_category(cat_df, asm_acc):
    category = cat_df[cat_df['asm_acc'] == asm_acc]['category'].values[0]
    if pd.isnull(category):
        category = 'NA'
    # end if
    return category
# end def


# == Proceed ==

# Read taxonomy file
tax_df = pd.read_csv(
    tax_fpath,
    sep='\t',
    dtype={
        'asm_acc': str,
        'taxid': pd.Int32Dtype(),
        'organism_name': str,
        'Species': str,
        'Genus': str,
        'Family': str,
        'Order': str,
        'Class': str,
        'Phylum': str,
        'Kingdom': str,
        'Domain': str,
    }
)

# Read categories file
cat_df = pd.read_csv(
    cat_fpath,
    sep='\t',
    dtype={
        'asm_acc': str,
        'category': pd.Int8Dtype(),
    }
)

# Count sequences
n_seqs = len(tuple(SeqIO.parse(in_fasta_fpath, 'fasta')))

# In order to select sequences quickly
tax_df.index = tax_df['asm_acc']

step = 1000
next_report_i = step
print(f'0/{n_seqs}', end='')
sys.stdout.flush()

cache_asm_acc = ''


with open(outfpath, 'wt') as outfile:

    seq_records = SeqIO.parse(in_fasta_fpath, 'fasta')

    i = 0
    for i, seq_record in enumerate(seq_records):
        asm_acc = parse_asm_acc(seq_record.id)

        # No need to create `taxonomy_str` and `category` from scratch each time
        if asm_acc != cache_asm_acc:
            taxonomy_str = get_taxonomy(asm_acc, tax_df)
            category = get_category(cat_df, asm_acc)
            cache_asm_acc = asm_acc
        # end if

        # Form header for output fasta file
        seq_record.description = '{} {}{}{} category:{}'.format(
            seq_record.id,
            TAX_SEP, taxonomy_str, TAX_SEP,
            category
        )

        outfile.write(f'>{seq_record.description}\n{seq_record.seq}\n')

        if i + 1 == next_report_i:
            print(f'\r{i+1}/{n_seqs}', end=' '*10)
            sys.stdout.flush()
            next_report_i += step
        # end if
    # end for
    print(f'\r{i+1}/{n_seqs}', end=' '*10)
# end with

print('\nCompleted!')
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
