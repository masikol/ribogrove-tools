#!/usr/bin/env python3

# The script finds repeats in gene sequences using RepeatFinder
#   https://github.com/deprekate/RepeatFinder

## Command line arguments

### Input files:
# 1. `-f / --in-fasta-file` -- an input fasta file of gene sequences.
#   This file is the output of the script `drop_aberrant_genes.py`. Mandatory.

### Output files:
# 1. `-o / --outfile` -- an output TSV file.
# The output file contains the following columns:
# - `seqID` -- RybaSom sequence identifier;
# - `r1_start`, `r1_end`, `r2_start`, `r2_end` -- cordinates of repeats within RybaSom sequences;
# - `rep_len` -- repeat length;
# - `rep_seq` -- repeat sequence;

## Dependencies:
# 1. RepeatFinder must be installed. See https://github.com/deprekate/RepeatFinder

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import argparse
from typing import Tuple

import repeatfinder as rf # https://github.com/deprekate/RepeatFinder
from Bio import SeqIO


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--in-fasta-file',
    help='input fasta file',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='TSV file for saving information about discovered repeats in genes sequences',
    required=True
)

args = parser.parse_args()


# For convenience
seqs_fpath = os.path.abspath(args.in_fasta_file)
# conserved_regions_fpath = os.path.abspath(args.conserved_regions_fasta)
outfpath = os.path.abspath(args.outfile)


if not os.path.exists(seqs_fpath):
    print(f'Error: file `{seqs_fpath}` does not exist!')
    sys.exit(1)
# end if

# Create output directory if needed
if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

print(seqs_fpath)
print()


# Some values for status messages
next_report = 499
inc = 500

num_seqs = len(tuple(SeqIO.parse(seqs_fpath, 'fasta'))) # tally input sequences
seq_records = SeqIO.parse(seqs_fpath, 'fasta') # read input sequences


def get_repeat_len(repeat_out: Tuple[int, int, int, int]):
    # Function for calculating length of a repeat
    return repeat_out[1] - repeat_out[0] + 1
# end def get_repeat_len


# == Proceed ==

with open(outfpath, 'wt') as outfile:

    # Write header
    outfile.write('seqID\tgene_len\tr1_start\tr1_end\tr2_start\tr2_end\trep_len\trep_seq\n')


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

            rep_len = get_repeat_len(r) # get length of repeat

            # Write output line
            outfile.write(f'{record.id}\t{len(record.seq)}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{rep_len}\t{rep_seq}\n')
        # end for
    # end for
# end with

print(f'\r{i+1}/{num_seqs}')
print('Completed!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
