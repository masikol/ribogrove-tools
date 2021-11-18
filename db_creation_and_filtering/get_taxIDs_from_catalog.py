#!/usr/bin/env python3

# The script maps Assembly IDs to Taxonomy IDs using elink utility
#   (https://www.ncbi.nlm.nih.gov/books/NBK25497/).
# Requires Internet connection.

## Command line arguments
### Input files:
# 1. `-i / --assm-acc-file` -- a TSV file of 4 columns: (`ass_id`, `gi_number`, `acc`, `title`).
#   This file is the output of the script `merge_assID2acc_and_remove_WGS.py`. Mandatory.
# 2. `-f / --all-fasta-file` -- a fasta file with all extracted genes sequences.
#   This file is the output of the script `extract_16S.py`. Mandatory.

### Output files:
# 1. `--per-genome-outfile` -- an output TSV file mapping Assembly IDs to taxIDs.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import re
import sys
import time
import gzip
import argparse
import subprocess as sp
from typing import List, Dict

import pandas as pd
from Bio import Entrez
from Bio import SeqIO


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '-c',
    '--refseq-catalog-file',
    help='the RefSeq catalog file (gzipped or not)',
    required=True
)

# Output files

parser.add_argument(
    '--per-genome-outfile',
    help='output file mapping Assembly IDs to taxIDs',
    required=True
)

args = parser.parse_args()


# For convenience
assm_acc_fpath = os.path.abspath(args.assm_acc_file)
# fasta_seqs_fpath = os.path.abspath(args.all_fasta_file)
refseq_catalog_fpath = os.path.abspath(args.refseq_catalog_file)
per_genome_outfpath = os.path.abspath(args.per_genome_outfile)
# email = args.email


# Check existance of all input files and dependencies
for fpath in (assm_acc_fpath, refseq_catalog_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# end for

# Create output directories if needed
if not os.path.isdir(os.path.dirname(per_genome_outfpath)):
    try:
        os.makedirs(os.path.dirname(per_genome_outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(per_genome_outfpath)}`')
        sys.exit(1)
    # end try
# end if

print(assm_acc_fpath)
print(refseq_catalog_fpath)
print()


# Read Assembly IDs, ACCESSION.VERSION's dataframe
assm_acc_df = pd.read_csv(
    assm_acc_fpath,
    sep='\t'
)

# Create tuple of Assembly IDs
ass_ids = tuple(
    set(
        assm_acc_df['ass_id']
    )
)

# Create a set of ACCESSION.VERSIONs
all_accs_set = set(
    assm_acc_df['acc']
)


# Read the catalog file
if refseq_catalog_fpath.endswith('.gz'):
    open_func = gzip.open
else:
    open_func = open
# end if

print(f'Reading the catalog file `{refseq_catalog_fpath}`')
acc_2_taxID_dict = dict()
with open_func(refseq_catalog_fpath, 'rt') as catalog_file:
    taxID_column_index = 0
    acc_column_index = 2
    separator = '\t'

    for line in catalog_file:
        col_vals = line.split(separator)
        if col_vals[acc_column_index] in all_accs_set:
            acc_2_taxID_dict[col_vals[acc_column_index]] = col_vals[taxID_column_index]
        # end if
    # end for
# end with
print('done\n')

del all_accs_set # don't need it any more

status_inc = 100
next_status_done_num = status_inc


# == Proceed ==

print('Mapping Assembly IDs to Taxonomy IDs...')

with open(per_genome_outfpath, 'wt') as per_genome_outfile:

    # Write headers to output files
    per_genome_outfile.write('ass_id\taccs\ttaxID\n')

    # Iterate over Assembly IDs
    for i, ass_id in enumerate(ass_ids):

        if i + 1 == next_status_done_num:
            print(f'\r{i+1}/{len(ass_ids)} taxIDs mapped', end=' '*10)
            next_status_done_num += status_inc
        # end if

        # Get all ACCESSION.VERSION's of current assembly
        accs = tuple(
            assm_acc_df[assm_acc_df['ass_id'] == ass_id]['acc']
        )

        tax_ID_set = set()
        for acc in accs:
            tax_ID_set.add(acc_2_taxID_dict[acc])
        # end if

        # Grt first-best taxID
        tax_id = next(iter(tax_ID_set))

        # Report multiple taxIDs -- let a user to check by him(her)self
        if len(tax_ID_set) != 1:
            print(f'\nMultiple TaxIDs for Assembly ID {ass_id}: {", ".join(tax_ID_set)}')
            print(f'Using this taxID: {tax_id}')
            # print(f'{ass_id} ({accs}): {tax_ID_set}')
        # end if

        # Write to per-genome output file
        per_genome_outfile.write(f'{ass_id}\t{";".join(accs)}\t{tax_id}\n')
    # end for
# end with

print(f'\r{i+1}/{len(ass_ids)} taxIDs mapped', end='\n\n')

print('Completed!')
print(per_genome_outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
