#!/usr/bin/env python3

# TODO: add description

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import argparse
from typing import Tuple

import numpy as np
import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--in-short-out-tsv',
    help='input .short.out.tsv file',
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
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
