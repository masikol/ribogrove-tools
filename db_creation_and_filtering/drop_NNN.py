#!/usr/bin/env python3

# The script discards gene sequences, which originate from the genomes containing at least
#   3 Ns (undefined bases) in a row in theis sequences.

## Command line arguments

### Input files:
# 1. `-i / --assm-acc-file` -- a TSV file of 4 columns: (`ass_id`, `gi_number`, `acc`, `title`).
#   This file is the output of the script `merge_assID2acc_and_remove_WGS.py`. Mandatory.
# 2. `-f / --all-fasta-file` -- a fasta file with all extracted genes sequences.
#   This file is the output of the script `extract_16S.py`. Mandatory.
# 3. `-c / --categories-file` -- a TSV file with genome categories.
#   This file is the output of the script `assign_genome_categories.py`. Mandatory.

### Output files:
# 1. `--out-fasta-file` -- a fasta file containing no sequences with NNN. Mandatory.
# 2. `--out-stats-file` -- a TSV file containing per-replicon statistics for
#   the output fasta file ("no NNN" one). Mandatory.
# 3. `--NNN-outfile` -- a fasta file containing sequences with NNN. Mandatory.

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import argparse

import pandas as pd
from Bio import SeqIO

from gene_seqs_2_stats import gene_seqs_2_stats


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
    '--out-fasta-file',
    help='output fasta file containing genes sequences without NNNs',
    required=True
)

parser.add_argument(
    '--out-stats-file',
    help='output per-replicon statistics file of sequences without NNNs',
    required=True
)

parser.add_argument(
    '--NNN-outfile',
    help='output fasta file containing genes sequences with NNNs',
    required=True
)


args = parser.parse_args()


assm_acc_fpath = os.path.abspath(args.assm_acc_file)
seqs_fpath = os.path.abspath(args.all_fasta_file)
category_fpath = os.path.abspath(args.categories_file)
out_fasta_fpath = os.path.abspath(args.out_fasta_file)
out_stats_fpath = os.path.abspath(args.out_stats_file)
nnn_seqs_fpath = os.path.abspath(args.NNN_outfile)


# Check existance of all input files
for fpath in (assm_acc_fpath, seqs_fpath, category_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directories if needed
for some_dir in map(os.path.dirname, [out_fasta_fpath, nnn_seqs_fpath, out_stats_fpath]):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
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


with open(out_fasta_fpath, 'wt') as out_fasta_file, \
     open(nnn_seqs_fpath, 'wt') as outfile_nnn:

    # Iterate over genes sequences
    for i, record in enumerate(seq_records):

        if i == next_report:
            # Proint status message
            print(f'\r{i+1}/{num_seqs}', end=' ')
            next_report += inc
        # end if

        if not record.id in NNN_seqIDs:
            # Write no "no NNN" file
            out_fasta_file.write(f'>{record.description}\n{str(record.seq)}\n')
        else:
            # Write no "NNN" file
            nnn_count += 1
            outfile_nnn.write(f'>{record.description}\n{str(record.seq)}\n')
        # end if
    # end for
# end with

print(f'\r{i+1}/{num_seqs}\n')
print(f'{nnn_count} genes from genomes with NNN found')
print(out_fasta_fpath)
print(nnn_seqs_fpath)

# Calculate statistics
print('Calculating statistics')
gene_seqs_2_stats(out_fasta_fpath, assm_acc_fpath, out_stats_fpath)
print()
print(out_stats_fpath)

print('Completed!')
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
