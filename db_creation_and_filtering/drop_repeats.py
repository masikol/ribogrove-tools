#!/usr/bin/env python3

# Script removes gene sequences, which contain large repeats.
# Threshold (repeat length) for distinguishing "large" and "not large" repeats
#   is automatically defined.

# Input files:
# 1. `-i/--assm-acc-file` is output of script merge_assID2acc_and_remove_WGS.py.
#   It has 4 columns: ass_id, refseq_id, acc, title. `refseq_id` is GI number.
# 2. Fasta file with all extracted genes sequences (-f/--all-fasta-file).

# Output files:
# 1. Fasta file containing no sequences with NN (--out-fasta-file).
# 2. `--out-stats-file` is a TSV file containing per-replicon statisticsw for `--out-fasta-file`.
# 3. Fasta file containing sequences with NN (--NN-outfile).

import os
import sys
import argparse
from typing import Sequence, Dict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gene_seqs_2_stats import gene_seqs_2_stats


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--input-fasta-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-i',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '-r',
    '--repeats-file',
    help="""TSV file (with header) with
  coordinates and lengths of each repeat found by script find_repeats.py""",
    required=True
)

parser.add_argument(
    '-e',
    '--exception-seqIDs',
    help="""Text file with seqIDs (one per line) of sequences, which are exceptions:
exceptions contain repeats, but they should not be removed""",
    required=True
)


# Output files

parser.add_argument(
    '--out-fasta-file',
    help='output fasta file containing sequences of genes without repeats',
    required=True
)

parser.add_argument(
    '--seqs-with-repeats',
    help='output fasta file containing sequences of genes with repeats',
    required=True
)

parser.add_argument(
    '--out-stats-file',
    help='output per-replicon statistics file',
    required=True
)

# Params

parser.add_argument(
    '--repeat-len-threshold',
    help='repeat length threshold (int > 0). Sequences with repeats longer than this threshold will be removed',
    required=True
)


args = parser.parse_args()

# For convenience
in_fasta_fpath = os.path.abspath(args.input_fasta_file)
assm_acc_fpath = os.path.abspath(args.assm_acc_file)
repeats_fpath = os.path.abspath(args.repeats_file)
exception_seqIDs_fpath = ps.path.abspath(args.exception_seqIDs)
output_genes_fpath = os.path.abspath(args.out_fasta_file)
seqs_with_repeats_fpath = os.path.abspath(args.seqs_with_repeats)
output_genes_stats_fpath = os.path.abspath(args.out_stats_file)


# Validate repeat_len_threshold
try:
    repeat_len_threshold = int(args.repeat_len_threshold)
    if repeat_len_threshold < 0:
        raise ValueError
    # end if
except ValueError:
    print(f'Invalid value of --repeat-len-threshold: `{args.repeat_len_threshold}`')
    print('It must be integer > 0.')
    sys.exit(1)
# end try


# Check existance of all input files
for fpath in (in_fasta_fpath, assm_acc_fpath, repeats_fpath, exception_seqIDs_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directories if needed
for some_dir in map(os.path.dirname, [output_genes_fpath, seqs_with_repeats_fpath, output_genes_stats_fpath]):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end if

exception_seqIDs = set(
    map(
        str.strip,
        open(exception_seqIDs_fpath, 'rt').readlines()
    )
)

print(in_fasta_fpath)
print(assm_acc_fpath)
print(repeats_fpath)
print(f'Repeats length threshold = {repeat_len_threshold}')
print()


# == Proceed ==

# Read sequence records
seq_records = tuple(SeqIO.parse(in_fasta_fpath, 'fasta'))

# Read repets dataframe
repeats_df = pd.read_csv(repeats_fpath, sep='\t')

# Select seqIDs with long repeats.
# And "substact" exception from them.
seqIDs_with_large_repeats = set(
    repeats_df[repeats_df['rep_len'] > repeat_len_threshold]['seqID']
) - exception_seqIDs

print(f'Found {len(seqIDs_with_large_repeats)} sequences with repeats longer than {repeat_len_threshold} bp')
print(f'{len(seq_records) - len(seqIDs_with_large_repeats)} sequences remain')
print(f'Writing sequences without repeats to file `{output_genes_fpath}`')
print(f'Writing sequences with repeats to file `{seqs_with_repeats_fpath}`')


# Write sequences to both files
with open(output_genes_fpath, 'wt') as output_genes_file, open(seqs_with_repeats_fpath, 'wt') as seqs_with_repeats_file:
    for seq_record in seq_records:
        if not seq_record.id in seqIDs_with_large_repeats:
            output_genes_file.write(f'>{seq_record.description}\n{seq_record.seq}\n')
        else:
            seqs_with_repeats_file.write(f'>{seq_record.description}\n{seq_record.seq}\n')
        # end if
    # end for
# end with

print('Done')
print(output_genes_fpath)
print(seqs_with_repeats_fpath)


print()

# Calculate per-replicon statistics
print('Calculating per-replicon statistics for remaining sequences')
gene_seqs_2_stats(output_genes_fpath, assm_acc_fpath, output_genes_stats_fpath)
print()

print('\nCompleted!')
print(output_genes_stats_fpath)
