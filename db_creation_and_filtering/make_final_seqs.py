#!/usr/bin/env python3

# This script produces a final (yet unannotated) file of target genes sequences.
# It takes "all collected" file and several filter files (of seqIDs, one eper line)
#   and filters these seqIDs out.
# Moreover, it takes a blacklist and a whitelist file (see `ad_hoc` directory).
#   Again, seqIDs, one per line.
# Whitelist: sequences, which the pipeline discards but you are sure they 
#   should not be discarded.
# Blacklist: sequences, which the pipeline doesn't discard but you are sure they 
#   should be discarded.

## Command line arguments
### Input files:
# 1. `-i / --all-seqs-file` -- an input fasta file of gene sequences.
#   This file is the output of the script `extract_16S.py`. Mandatory.
# 2. `--NNN-fail-seqIDs` -- a file of seqIDs which didn't pass NNN filter, one per line.
#   This file is the output of the script `find_NNN.py`. Mandatory.
# 3. `--aberrant-seqIDs` -- a file of seqIDs which didn't pass "find aberrant genes" filter, one per line.
#   This file is the output of the script `find_aberrant_genes.py`. Mandatory.
# 4. `--repeats-fail-seqIDs` -- a file of seqIDs which didn't pass "repeats" filter, one per line.
#   This file is the output of the script `find_repeats.py`. Mandatory.
# 5. `--blacklist-seqIDs` -- a file of seqIDs you know should be discarded. Mandatory.
# 6. `--whitelist-seqIDs` -- a file of seqIDs you know should not be discarded. Mandatory.
# 7. `-a / --assm-acc-file` -- a TSV file of 4 columns: (`ass_id`, `gi_number`, `acc`, `title`).
#   This file is the output of the script `merge_assIDs_and_accs.py`. Mandatory.

### Output files:
# 1. `--out-fasta-file` -- output fasta file of final sequences.
#   Mandatory.
# 2. `--out-stats-file` -- a file with per-replicon statistics of filtered 16S genes:
#   how many genes, minimum/maximum length etc.
#   Mandatory.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import argparse

import pandas as pd
from Bio import SeqIO

from gene_seqs_2_stats import gene_seqs_2_stats
from read_and_filter_fasta import read_and_filter_fasta


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--all-seqs-file',
    help='input fasta file of all collected sequences',
    required=True
)

parser.add_argument(
    '--NNN-fail-seqIDs',
    help='input file of seqIDs which didn\'t pass NNN filter',
    required=True
)

parser.add_argument(
    '--aberrant-seqIDs',
    help='input file of seqIDs which didn\'t pass the "find aberrant" filter',
    required=True
)

parser.add_argument(
    '--repeats-fail-seqIDs',
    help='input file of seqIDs which didn\'t pass the repeat filter',
    required=True
)

parser.add_argument(
    '--blacklist-seqIDs',
    help='a file of seqIDs you know should be discarded',
    required=True
)

parser.add_argument(
    '--whitelist-seqIDs',
    help='a file of seqIDs you know should not be discarded',
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
    help='output fasta file of final sequences',
    required=True
)

parser.add_argument(
    '--out-stats-file',
    help='a file with per-replicon statistics of filtered 16S genes',
    required=True
)

args = parser.parse_args()


# For convenience
input_seqs_fpath = os.path.abspath(args.all_seqs_file)
NNN_fail_fpath = os.path.abspath(args.NNN_fail_seqIDs)
aberrant_fpath = os.path.abspath(args.aberrant_seqIDs)
repeats_fail_fpath = os.path.abspath(args.repeats_fail_seqIDs)
blacklist_fpath = os.path.abspath(args.blacklist_seqIDs)
whitelist_fpath = os.path.abspath(args.whitelist_seqIDs)
assm_acc_fpath = os.path.abspath(args.assm_acc_file)
out_fasta_fpath = os.path.abspath(args.out_fasta_file)
out_stats_fpath = os.path.abspath(args.out_stats_file)


fpaths_to_check = (
    input_seqs_fpath,
    NNN_fail_fpath,
    aberrant_fpath,
    repeats_fail_fpath,
    assm_acc_fpath,
    blacklist_fpath,
    whitelist_fpath,
)
for f in fpaths_to_check:
    if not os.path.exists(f):
        print(f'Error: file `{f}` does not exist!')
        sys.exit(1)
    # end if
# end for
del f, fpaths_to_check

# Create output directory if needed
for d in map(os.path.dirname, [out_fasta_fpath, out_stats_fpath]):
    if not os.path.isdir(d):
        try:
            os.makedirs(d)
        except OSError as err:
            print(f'Error: cannot create directory `{d}`')
            sys.exit(1)
        # end try
    # end if
# end for

print(input_seqs_fpath)
print(NNN_fail_fpath)
print(aberrant_fpath)
print(repeats_fail_fpath)
print(assm_acc_fpath)
print(blacklist_fpath)
print(whitelist_fpath)
print()


print("Started")
print('1. Filtering sequences...')

blacklist = set(
    pd.read_csv(blacklist_fpath, sep='\t')['seqID']
)
whitelist = set(
    pd.read_csv(whitelist_fpath, sep='\t')['seqID']
)

# Read input sequences
final_seq_records = read_and_filter_fasta(
    input_seqs_fpath,
    filter_fpaths=[
        NNN_fail_fpath,
        aberrant_fpath,
        repeats_fail_fpath,
    ],
    blacklist=blacklist,
    whitelist=whitelist
)

with open(out_fasta_fpath, 'wt') as out_fasta_file:
    SeqIO.write(final_seq_records, out_fasta_file, 'fasta')
# end with
print(out_fasta_fpath)
print('Done')


print()

# Calculate per-replicon statistics
print('2. Calculating per-replicon statistics for filtered sequences')
gene_seqs_2_stats(out_fasta_fpath, assm_acc_fpath, out_stats_fpath)
print()
print('Done')
print(out_stats_fpath)

print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
