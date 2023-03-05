#!/usr/bin/env python3

# The script records which gene sequences didn’t pass the `ribotyper` filter,
#   i.e. sequences with the following unexpected features:
#   "NoHits", "UnacceptableModel", "MinusStrand", "LowScore", "LowCoverage".

## Command line arguments

### Input files:
# 1. `-i / --in-short-out-tsv` -- input .short.out.tsv file
#   outputted by the script `check_seqs_with_ribotyper.py`.
#   Mandatory.

### Output files:
# 1. `-o / --out-fail-file` -- output file of seqIDs which don't pass the filter,
#   one per line.
#   Mandatory


import os
from src.rg_tools_time import get_time

print(
    '\n|=== {} STARTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)


# == Parse arguments ==
import argparse

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--in-short-out-tsv',
    help="""input .short.out.tsv file
    outputted by the script `check_seqs_with_ribotyper.py`""",
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--out-fail-file',
    help='output file of seqIDs which don\'t pass the filter, one per line',
    required=True
)

args = parser.parse_args()


# == Import them now ==
import sys
from typing import Tuple

import numpy as np
import pandas as pd


# For convenience
in_short_out_fpath = os.path.abspath(args.in_short_out_tsv)
out_fail_fpath = os.path.abspath(args.out_fail_file)


if not os.path.exists(in_short_out_fpath):
    print(f'Error: file `{in_short_out_fpath}` does not exist!')
    sys.exit(1)
# end if

# Create output directory if needed
if not os.path.isdir(os.path.dirname(out_fail_fpath)):
    try:
        os.makedirs(os.path.dirname(out_fail_fpath))
    except OSError as err:
        print(
            'Error: cannot create directory `{}`' \
                .format(os.path.dirname(out_fail_fpath))
        )
        sys.exit(1)
    # end try
# end if


print(in_short_out_fpath)
print()


FAIL_FEATURES = (
    'NoHits',
    'UnacceptableModel',
    'MinusStrand',
    'LowScore',
    'LowCoverage',
)


def find_failed_seqIDs(in_short_out_fpath):
    short_out_df = pd.read_csv(in_short_out_fpath, sep='\t')
    short_out_df['custom_fail'] = np.repeat(False, short_out_df.shape[0])
    short_out_df = short_out_df.apply(set_custom_fail, axis=1)
    return tuple(
        short_out_df[short_out_df['custom_fail'] == True]['target']
    )
# end def

def set_custom_fail(row):
    global FAIL_FEATURES

    row_features = row['unexpected_features']
    if row_features == '-':
        return row
    # end if

    for fail_feature in FAIL_FEATURES:
        if fail_feature in row_features:
            row['custom_fail'] = True
            return row
        # end if
    # end for
    return row
# end def


# == Proceed ==

print(
    'Searching sequences having following features: {}'.format(FAIL_FEATURES)
)

failed_seqIDs = find_failed_seqIDs(in_short_out_fpath)

with open(out_fail_fpath, 'w') as outfile:
    outfile.write('\n'.join(failed_seqIDs))
    outfile.write('\n')
# end with

print('Completed!')
print(out_fail_fpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
