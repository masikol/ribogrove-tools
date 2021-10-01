#!/usr/bin/env python3

import os
import sys
import argparse

import pandas as pd
from Bio import SeqIO

from gene_seqs_2_stats import gene_seqs_2_stats


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--in-fasta-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-t',
    '--tblout',
    help='output .tblout file of cmscan',
    required=True
)

parser.add_argument(
    '-a',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)


# Output files

parser.add_argument(
    '--out-fasta-file',
    help='output fasta file',
    required=True
)

parser.add_argument(
    '--trunc-fasta-file',
    help='output fasta file for truncated sequences',
    required=True
)

parser.add_argument(
    '--out-stats-file',
    help='output per-replicon statistics file',
    required=True
)



args = parser.parse_args()

# For convenience
fasta_seqs_fpath = os.path.abspath(args.in_fasta_file)
tblout_fpath = os.path.abspath(args.tblout)
assm_acc_fpath = os.path.abspath(args.assm_acc_file)
out_fasta_fpath = os.path.abspath(args.out_fasta_file)
trunc_fasta_fpath = os.path.abspath(args.trunc_fasta_file)
out_stats_fpath = os.path.abspath(args.out_stats_file)


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, tblout_fpath, assm_acc_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for


# Create output directories if needed
for some_dir in map(os.path.dirname, [out_fasta_fpath, out_stats_fpath]):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end for

print(fasta_seqs_fpath)
print(trunc_fasta_fpath)
print(tblout_fpath)
print(assm_acc_fpath)
print()

tblout_df = pd.read_csv(
    tblout_fpath,
    sep='\t'
)


input_seq_records = tuple(SeqIO.parse(fasta_seqs_fpath, 'fasta'))

seqIDs_to_keep = set(tblout_df[tblout_df['trunc'] == 'no']['query_name'])

print(f'{len(input_seq_records) - len(seqIDs_to_keep)} sequences are truncated')
print(f'{len(seqIDs_to_keep)} sequences remain')
print(f'Writing full-length sequences to file `{out_fasta_fpath}`...')
print(f'Writing truncated sequences to trash file `{trunc_fasta_fpath}`...')


with open(out_fasta_fpath, 'wt') as out_fasta_file, open(trunc_fasta_fpath, 'wt') as trunc_fasta_file:
    for seq_record in input_seq_records:
        if seq_record.id in seqIDs_to_keep:
            out_fasta_file.write(f'>{seq_record.description}\n{seq_record.seq}\n')
        else:
            trunc_fasta_file.write(f'>{seq_record.description}\n{seq_record.seq}\n')
        # end if
    # end for
# end with
print('Done')
print(out_fasta_fpath)
print(trunc_fasta_fpath)


# Calculate statistics
print('Calculating statistics for result full-length sequences')
gene_seqs_2_stats(out_fasta_fpath, assm_acc_fpath, out_stats_fpath)
print()
print(out_stats_fpath)

print('Completed!')
