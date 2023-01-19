#!/usr/bin/env python3

# The script finds gene sequences, which originate from the genomes containing at least
#   3 Ns (undefined bases) in a row in theis sequences.

## Command line arguments

### Input files:
# 1. `-i / --assm-acc-file` -- a TSV file of 4 columns: (`ass_id`, `gi_number`, `acc`, `title`).
#   This file is the output of the script `merge_assIDs_and_accs.py`. Mandatory.
# 2. `-f / --all-fasta-file` -- a fasta file with all collected genes sequences.
#   This file is the output of the script `extract_16S.py`. Mandatory.
# 3. `-c / --categories-file` -- a TSV file with genome categories.
#   This file is the output of the script `assign_genome_categories.py`. Mandatory.

### Output files:
# 1. `--out-fail-file` -- output file of seqIDs which don't pass the filter, one per line.
#   Mandatory.

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
    '-i',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '-f',
    '--all-fasta-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-c',
    '--categories-file',
    help='TSV file of genome categories info',
    required=True
)

# Output files

parser.add_argument(
    '--out-fail-file',
    help='output file of seqIDs which don\'t pass the filter, one per line',
    required=True
)


args = parser.parse_args()


assm_acc_fpath = os.path.abspath(args.assm_acc_file)
seqs_fpath = os.path.abspath(args.all_fasta_file)
category_fpath = os.path.abspath(args.categories_file)
out_fail_fpath = os.path.abspath(args.out_fail_file)


# Check existance of all input files
for fpath in (assm_acc_fpath, seqs_fpath, category_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directories if needed
if not os.path.isdir(os.path.dirname(out_fail_fpath)):
    try:
        os.makedirs(os.path.dirname(out_fail_fpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(out_fail_fpath)}`')
        sys.exit(1)
    # end try
# end if


nnn_count = 0

next_report = 499
inc = 500


num_seqs = len(tuple(SeqIO.parse(seqs_fpath, 'fasta'))) # tally sequences
seq_records = SeqIO.parse(seqs_fpath, 'fasta') # read sequences


# == Proceed ==

category_df = pd.read_csv(
    category_fpath,
    sep='\t'
)

NNN_seqIDs = set(
    category_df[category_df['contains_NNN'] == 1]['seqID']
)


with open(out_fail_fpath, 'wt') as out_fail_file:
    # Iterate over genes sequences
    for i, record in enumerate(seq_records):

        if i == next_report:
            # Proint status message
            print(f'\r{i+1}/{num_seqs}', end=' ')
            next_report += inc
        # end if

        if record.id in NNN_seqIDs:
            nnn_count += 1
            out_fail_file.write(f'{record.id}\n')
        # end if
    # end for
# end with

print(f'\r{i+1}/{num_seqs}\n')
print(f'{nnn_count} genes from genomes with NNN found')
print(out_fail_fpath)

print('Completed!')
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
