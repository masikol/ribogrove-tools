#!/usr/bin/env python3

import os
import sys
import argparse

import polars as pl

from src.rg_tools_time import get_time


sys.stderr.write(
    '\n|=== {} STARTING SCRIPT `{}` ===|\n\n'.format(
        get_time(), os.path.basename(__file__)
    )
)


parser = argparse.ArgumentParser()
parser.add_argument('-c', '--current-hash-table', required=True)
parser.add_argument('-p', '--prev-hash-table', required=True)
parser.add_argument('-o', '--out', required=True)
args = parser.parse_args()


# == Import them now ==

curr_df_fpath = os.path.abspath(args.current_hash_table)
prev_df_fpath = os.path.abspath(args.prev_hash_table)
out_fpath = os.path.abspath(args.out)

for fpath in (curr_df_fpath, prev_df_fpath):
    if not os.path.exists(fpath):
        sys.stderr.write('Error: input file `{}` does not exist!\n'.format(fpath))
        sys.exit(1)
    # end if
# end for

out_dir = os.path.dirname(out_fpath)
if out_dir and not os.path.isdir(out_dir):
    try:
        os.makedirs(out_dir)
    except OSError as err:
        sys.stderr.write('Error: cannot create directory `{}`\n'.format(out_dir))
        sys.exit(1)
    # end try
# end if


# >>> Proceed >>>

curr_df = pl.read_csv(curr_df_fpath, separator='\t') \
    .rename({'asm_seqID_hash': 'curr_hash'})
prev_df = pl.read_csv(prev_df_fpath, separator='\t') \
    .rename({'asm_seqID_hash': 'prev_hash'})

diff_df = curr_df.join(prev_df, on='asm_acc', how='left')
diff_df = diff_df.filter(
    pl.col('curr_hash') != pl.col('prev_hash')
)

diff_df.write_csv(out_fpath, separator='\t')



sys.stderr.write(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n\n'.format(
        get_time(), os.path.basename(__file__)
    )
)

sys.exit(0)
