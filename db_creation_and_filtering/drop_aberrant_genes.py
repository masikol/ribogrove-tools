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
    '--input-fasta-file',
    help='fasta file of SSU gene sequences',
    required=True
)

# parser.add_argument(
#     '-t',
#     '--tblout',
#     help='output .tblout file of cmscan',
#     required=True
# )

parser.add_argument(
    '-a',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '--non-aberrant-seqIDs',
    help='txt file with seqIDs (one per line) of non-aberrant sequences',
    required=True
)

parser.add_argument(
    '--aberrant-seqIDs',
    help='txt file with seqIDs (one per line) of aberrant sequences',
    required=True
)



# Output files

parser.add_argument(
    '--non-aberrant-fasta-file',
    help='output fasta file',
    required=True
)

parser.add_argument(
    '--aberrant-fasta-file',
    help='output fasta file for aberrant sequences',
    required=True
)

parser.add_argument(
    '--out-stats-file',
    help='output per-replicon statistics file',
    required=True
)



args = parser.parse_args()

# For convenience
fasta_seqs_fpath = os.path.abspath(args.input_fasta_file)
# tblout_fpath = os.path.abspath(args.tblout)
assm_acc_fpath = os.path.abspath(args.assm_acc_file)
non_aberrant_seqIDs_fpath = os.path.abspath(args.non_aberrant_seqIDs)
aberrant_seqIDs_fpath = os.path.abspath(args.aberrant_seqIDs)
non_aberrant_fasta_fpath = os.path.abspath(args.non_aberrant_fasta_file)
aberrant_fasta_fpath = os.path.abspath(args.aberrant_fasta_file)
out_stats_fpath = os.path.abspath(args.out_stats_file)


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, assm_acc_fpath, non_aberrant_seqIDs_fpath, aberrant_seqIDs_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for


# Create output directories if needed
for some_dir in map(os.path.dirname, [non_aberrant_fasta_fpath, out_stats_fpath, aberrant_fasta_fpath]):
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
print(non_aberrant_seqIDs_fpath)
print(aberrant_seqIDs_fpath)
# print(tblout_fpath)
print(assm_acc_fpath)
print()

# tblout_df = pd.read_csv(
#     tblout_fpath,
#     sep='\t'
# )


input_seq_records = tuple(SeqIO.parse(fasta_seqs_fpath, 'fasta'))

non_aberrant_seqIDs = set(
    map(
        str.strip,
        open(non_aberrant_seqIDs_fpath, 'rt').readlines()
    )
)

aberrant_seqIDs = set(
    map(
        str.strip,
        open(aberrant_seqIDs_fpath, 'rt').readlines()
    )
)

# seqIDs_to_keep = set(tblout_df[tblout_df['trunc'] == 'no']['query_name'])

print(f'{len(aberrant_seqIDs)} sequences are considered as aberrant')
print(f'{len(non_aberrant_seqIDs)} sequences remain non-aberrant')
print(f'Writing non-aberrant sequences to file `{non_aberrant_fasta_fpath}`...')
print(f'Writing aberrant sequences to trash file `{aberrant_fasta_fpath}`...')


with open(non_aberrant_fasta_fpath, 'wt') as non_aberrant_fasta_file, open(aberrant_fasta_fpath, 'wt') as aberrant_fasta_file:
    for seq_record in input_seq_records:
        if seq_record.id in non_aberrant_seqIDs:
            non_aberrant_fasta_file.write(f'>{seq_record.description}\n{seq_record.seq}\n')
        elif seq_record.id in aberrant_seqIDs:
            aberrant_fasta_file.write(f'>{seq_record.description}\n{seq_record.seq}\n')
        else:
            print(f"""\nWarning: seqID `{seq_record.id}` is present neither in non-aberrant seqIDs
  nor in aberrant seqIDs. Please, check your input fasta file and files
  `--non-aberrant-seqIDs` and `--aberrant-seqIDs`""")
        # end if
    # end for
# end with
print('Done')
print(non_aberrant_fasta_fpath)
print(aberrant_fasta_fpath)


# Calculate statistics
print(f'Calculating statistics for result non-aberrant sequences in file `{non_aberrant_fasta_fpath}`')
gene_seqs_2_stats(non_aberrant_fasta_fpath, assm_acc_fpath, out_stats_fpath)
print()
print(out_stats_fpath)

print('Completed!')
