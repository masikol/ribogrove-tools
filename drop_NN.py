#!/usr/bin/env python3

# Script removes gene sequences, which containt at least 2 N's in theirs sequences

# Input files:
# 1. `-i/--assm-acc-file` is output of script merge_assID2acc_and_remove_WGS.py.
#   It has 4 columns: ass_id, refseq_id, acc, title. `refseq_id` is GI number.
# 2. Fasta file with all extracted genes sequences (-f/--all-fasta-file).

# Output files:
# 1. Fasta file containing no sequences with NN (--out-fasta-file).
# 2. `--out-stats-file` is a TSV file containing per-replicon statisticsw for `--out-fasta-file`.
# 3. Fasta file containing sequences with NN (--NN-outfile).


import os
import re
import sys
import argparse

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

# Output files

parser.add_argument(
    '--out-fasta-file',
    help='output fasta file containing genes sequences without NN\'s',
    required=True
)

parser.add_argument(
    '--out-stats-file',
    help='output per-replicon statistics file',
    required=True
)

parser.add_argument(
    '--NN-outfile',
    help='output fasta file containing genes sequences with NN\'s',
    required=True
)


args = parser.parse_args()


assm_acc_fpath = os.path.abspath(args.assm_acc_file)
seqs_fpath = os.path.abspath(args.all_fasta_file)
out_fasta_fpath = os.path.abspath(args.out_fasta_file)
out_stats_fpath = os.path.abspath(args.out_stats_file)
nn_seqs_fpath = os.path.abspath(args.NN_outfile)


# Check existance of all input files and dependencies
for fpath in (assm_acc_fpath, seqs_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directories if needed
for some_dir in map(os.path.dirname, [out_fasta_fpath, nn_seqs_fpath, out_stats_fpath]):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end if


nn_pattern = r'NN' # pattern for searching NN
nn_count = 0

next_report = 499
inc = 500


num_seqs = len(tuple(SeqIO.parse(seqs_fpath, 'fasta'))) # tally sequences
seq_records = SeqIO.parse(seqs_fpath, 'fasta') # read sequences


# == Proceed ==

with open(out_fasta_fpath, 'wt') as out_fasta_file, \
     open(nn_seqs_fpath, 'wt') as outfile_nn:

    # Iterate over genes sequences
    for i, record in enumerate(seq_records):

        if i == next_report:
            # Proint status message
            print(f'\r{i+1}/{num_seqs}', end=' ')
            next_report += inc
        # end if

        if re.search(nn_pattern, str(record.seq)) is None:
            # Write no "no NN" file
            out_fasta_file.write(f'>{record.description}\n{str(record.seq)}\n')
        else:
            # Write no "NN" file
            nn_count += 1
            outfile_nn.write(f'>{record.description}\n{str(record.seq)}\n')
        # end if
    # end for
# end with

print(f'\r{i+1}/{num_seqs}\n')
print(f'{nn_count} NN-sequences found')
print(out_fasta_fpath)
print(nn_seqs_fpath)

# Calculate statistics
print('Calculating statistics')
gene_seqs_2_stats(out_fasta_fpath, assm_acc_fpath, out_stats_fpath)
print()
print(out_stats_fpath)

print('Completed!')
