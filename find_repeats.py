#!/usr/bin/env python3

# Script finds repeats in SSU genes sequences using RepeatFinder.

# Input files:
# 1. Fasta file of genes sequences containing no NN (-f/--no-NN-fasta-file).
# 2. Fasta file of NR conserved regions from work
#    "How conserved are the conserved 16S-rRNA regions?"
#    (table 5, https://peerj.com/articles/3036/)
#    -c/--conserved-regions-fasta

# Output files:
# 1. TSV file (-o/--outfile) of following columns:
#   seqID;
#   r1_start, r1_end, r2_start, r2_end -- repeats' cordinates;
#   rep_len -- length of a repeat;
#   rep_seq -- sequence of a repeat;
#   conserv_<REGION_ID> -- columns indicating if repeat contains corresponding
#     conserved region from file -c/--conserved-regions-fasta.
#     1 (one) if contains, otherwise 0 (zero).

# Dependencies:
# 1. RepeatFinder must be installed. See https://github.com/deprekate/RepeatFinder


import os
import re
import sys
import argparse
from typing import Tuple

import repeatfinder as rf # https://github.com/deprekate/RepeatFinder
from Bio import SeqIO
from Bio import SeqUtils


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--no-NN-fasta-file',
    help='fasta file of SSU gene sequences containing no NN',
    required=True
)

parser.add_argument(
    '-c',
    '--conserved-regions-fasta',
    help="""fasta file of NR conserved regions from work
    "How conserved are the conserved 16S-rRNA regions?" (table 5, https://peerj.com/articles/3036/)""",
    required=True
)

parser.add_argument(
    '-o',
    '--outfile',
    help='TSV file for saving information about discovered repeats in genes sequences',
    required=True
)

args = parser.parse_args()


# For convenience
seqs_fpath = os.path.abspath(args.no_NN_fasta_file)
conserved_regions_fpath = os.path.abspath(args.conserved_regions_fasta)
outfpath = os.path.abspath(args.outfile)

# Check existance of all input files
for fpath in (seqs_fpath, conserved_regions_fpath):
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


conserved_seq_records = tuple(SeqIO.parse(conserved_regions_fpath, 'fasta'))

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
    outfile.write('seqID\tr1_start\tr1_end\tr2_start\tr2_end\trep_len\trep_seq\t')
    outfile.write('\t'.join(
        [f'conserv_{r.id}' for r in conserved_seq_records]
        ) + '\n'
    )

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

            # We assume that repeat contains no conserved regions
            conserv_list = ['0'] * len(conserved_seq_records)

            # Get repeat sequence
            rep_seq = str(record.seq)[r[0]-1 : r[1]]

            # Check if repeat contains conserved regions
            for i, conserv_record in enumerate(conserved_seq_records):
                search_list = SeqUtils.nt_search(rep_seq, str(conserv_record.seq))
                if len(search_list) > 1:
                    conserv_list[i] = '1'
                # end if
            # end for

            rep_len = get_repeat_len(r) # get length of repeat

            # Write output line
            outfile.write(f'{record.id}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{rep_len}\t{rep_seq}\t')
            outfile.write('{}\n'.format('\t'.join(conserv_list)))
        # end for
    # end for
# end with

print(f'\r{i+1}/{num_seqs}')
print('Completed!')
print(outfpath)
