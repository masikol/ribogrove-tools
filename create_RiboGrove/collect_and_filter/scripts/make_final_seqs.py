#!/usr/bin/env python3

# TODO: update description (--assm-acc-file)

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



import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import argparse

import pandas as pd
from Bio import SeqIO

import src.rg_tools_IO as rgIO


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
    '--ribotyper-fail-seqIDs',
    help='TODO: add help',
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
    '-m',
    '--replicon-map',
    help='TODO: add help',
    required=True
)

# Output files

parser.add_argument(
    '--out-fasta-file',
    help='output fasta file of final sequences',
    required=True
)

args = parser.parse_args()


# For convenience
input_seqs_fpath = os.path.abspath(args.all_seqs_file)
ribotyper_fail_fpath = os.path.abspath(args.ribotyper_fail_seqIDs)
aberrant_fpath = os.path.abspath(args.aberrant_seqIDs)
repeats_fail_fpath = os.path.abspath(args.repeats_fail_seqIDs)
blacklist_fpath = os.path.abspath(args.blacklist_seqIDs)
whitelist_fpath = os.path.abspath(args.whitelist_seqIDs)
replicon_map_fpath = os.path.abspath(args.replicon_map)
out_fasta_fpath = os.path.abspath(args.out_fasta_file)


fpaths_to_check = (
    input_seqs_fpath,
    ribotyper_fail_fpath,
    aberrant_fpath,
    repeats_fail_fpath,
    replicon_map_fpath,
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
if not os.path.isdir(os.path.dirname(out_fasta_fpath)):
    try:
        os.makedirs(os.path.dirname(out_fasta_fpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(out_fasta_fpath)}`')
        sys.exit(1)
    # end try
# end if


print(input_seqs_fpath)
print(ribotyper_fail_fpath)
print(aberrant_fpath)
print(repeats_fail_fpath)
print(replicon_map_fpath)
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
final_seq_records = rgIO.read_and_filter_fasta(
    input_seqs_fpath,
    filter_fpaths=[
        ribotyper_fail_fpath,
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


print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
