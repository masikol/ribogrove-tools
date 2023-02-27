#!/usr/bin/env python3

# The script finds repeats in gene sequences using RepeatFinder
#   https://github.com/deprekate/RepeatFinder

## Command line arguments

### Input files:
# 1. `-f / --in-fasta-file` -- an input fasta file of gene sequences.
#   This file is the output of the script `extract_16S.py`. Mandatory.
# 2. `--ribotyper-fail-seqIDs` -- a file of seqIDs which didn't pass ribotyper filter, one per line.
#   This file is the output of the script `find_ribotyper.py`. Mandatory.
# 3. `--aberrant-seqIDs` -- a file of seqIDs which didn't pass "find aberrant genes" filter, one per line.
#   This file is the output of the script `find_aberrant_genes.py`. Mandatory.

### Output files:
# 1. `--out-fail-file` -- a file of seqIDs which don't pass the filter, one per line.
#   Mandatory.
# 2. `--out-repeats-log` -- an output TSV listing all repets found in the input sequences.
# The file contains the following columns:
# - `seqID` -- RiboGrove sequence identifier;
# - `r1_start`, `r1_end`, `r2_start`, `r2_end` -- cordinates of repeats within RiboGrove sequences;
# - `rep_len` -- repeat length;
# - `rep_seq` -- repeat sequence;
#   Mandatory.

## Dependencies:
# 1. RepeatFinder must be installed. See https://github.com/deprekate/RepeatFinder

### Parameters:
# 1. `--repeat-len-threshold` -- a repeat length threshold.
#   Sequences having repeats longer than this value will be discarded. Mandatory.

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import argparse
from typing import Tuple

import repeatfinder as rf # https://github.com/deprekate/RepeatFinder
from Bio import SeqIO

import src.rg_tools_IO as rgIO


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--in-fasta-file',
    help='input fasta file',
    required=True
)

parser.add_argument(
    '--ribotyper-fail-seqIDs',
    help='input file of seqIDs which didn\'t pass ribotyper filter',
    required=True
)

parser.add_argument(
    '--aberrant-seqIDs',
    help='input file of seqIDs which didn\'t pass the "find aberrant" filter',
    required=True
)

# Output files

parser.add_argument(
    '--out-fail-file',
    help='output file of seqIDs which don\'t pass the filter, one per line',
    required=True
)

parser.add_argument(
    '--out-repeats-log',
    help='a file listing all repets found in the input sequences',
    required=True
)

# Params

parser.add_argument(
    '--repeat-len-threshold',
    help="""repeat length threshold (int > 0). Sequences with repeats longer than
  this threshold will be discarded""",
    required=True
)

args = parser.parse_args()


# For convenience
seqs_fpath = os.path.abspath(args.in_fasta_file)
ribotyper_fail_fpath = os.path.abspath(args.ribotyper_fail_seqIDs)
aberrant_fpath = os.path.abspath(args.aberrant_seqIDs)
out_fail_fpath = os.path.abspath(args.out_fail_file)
out_repeats_log_fpath = os.path.abspath(args.out_repeats_log)


for f in (seqs_fpath, ribotyper_fail_fpath, aberrant_fpath):
    if not os.path.exists(f):
        print(f'Error: file `{f}` does not exist!')
        sys.exit(1)
    # end if
# end for

# Create output directory if needed
for d in map(os.path.dirname, [out_fail_fpath, out_repeats_log_fpath]):
    if not os.path.isdir(d):
        try:
            os.makedirs(d)
        except OSError as err:
            print(f'Error: cannot create directory `{d}`')
            sys.exit(1)
        # end try
    # end if
# end for

# Parse and validate repeat_len_threshold
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

print(seqs_fpath)
print(ribotyper_fail_fpath)
print(aberrant_fpath)
print(f'Repeats length threshold = {repeat_len_threshold}')
print()


# Some values for status messages
next_report = 499
inc = 500

# Read input sequences
seq_records = rgIO.read_and_filter_fasta(
    seqs_fpath,
    filter_fpaths=[ribotyper_fail_fpath, aberrant_fpath,]
)

num_seqs = len(seq_records) # tally input sequences


def get_repeat_len(repeat_out: Tuple[int, int, int, int]):
    # Function for calculating length of a repeat
    return repeat_out[1] - repeat_out[0] + 1
# end def


# == Proceed ==

with open(out_repeats_log_fpath, 'wt') as out_log_file, \
     open(out_fail_fpath, 'wt')        as out_seqID_file:

    # Write header
    out_log_file.write('seqID\tgene_len\tr1_start\tr1_end\tr2_start\tr2_end\trep_len\trep_seq\n')


    # Iterate over input seq records
    for i, record in enumerate(seq_records):

        if i == next_report:
            # Primt status message
            print(f'\r{i+1}/{num_seqs}', end=' ')
            next_report += inc
        # end if

        # Find repeats
        repeats = rf.get_repeats(str(record.seq))

        # Record repeats
        for r in repeats:

            # Get repeat sequence
            rep_seq = str(record.seq)[r[0]-1 : r[1]]

            rep_len = get_repeat_len(r)

            # Write output line
            out_log_file.write(f'{record.id}\t{len(record.seq)}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{rep_len}\t{rep_seq}\n')

            if rep_len > repeat_len_threshold:
                out_seqID_file.write('{}\n'.format(record.id))
            # end if
        # end for
    # end for
# end with

print(f'\r{i+1}/{num_seqs}')
print('Completed!')
print(out_fail_fpath)
print(out_repeats_log_fpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
