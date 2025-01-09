#!/usr/bin/env python3

# The script processes extracted gene sequences with ribotyper program from Ribovore suite.

## Command line arguments

### Input files:
# 1. `-f / --in-fasta-file` -- input fasta file of 16S rRNA gene sequences.
#   Mandatory.

### Output files:
# 1. `-o / --outdir` -- output directory.
#   Mandatory.

### Dependencies:
# 1. `--ribotyper` -- path to ribotyper executable.
#   Mandatory.
# 2. `--acccept-file` -- file .accept for ribotyper.
#   It should be one of the followong files:
#   `data/archaea.accept`
#   `data/bacteria.accept`
#   Mandatory.

### "Cached" files:
# 1. `--prev-short-out-tsv` -- file `ribotyper_out.ribotyper.short.out` from the
#   previous RiboGrove release.
# 2. `--prev-long-out-tsv`` -- file `ribotyper_out.ribotyper.long.out` from the
#   previous RiboGrove release.

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
    '-f',
    '--in-fasta-file',
    help='input fasta file of SSU gene sequences',
    required=True
)

# "Cached" files

parser.add_argument(
    '--prev-short-out-tsv',
    help="""file `ribotyper_out.ribotyper.short.out` from the
    previous RiboGrove release""",
    required=False
)

parser.add_argument(
    '--prev-long-out-tsv',
    help="""file `ribotyper_out.ribotyper.long.out` from the
    previous RiboGrove release""",
    required=False
)

# Output files

parser.add_argument(
    '-o',
    '--outdir',
    help='output directory',
    required=True
)

# Dependencies

parser.add_argument(
    '--ribotyper',
    help='path to ribotyper executable',
    required=True
)

parser.add_argument(
    '--acccept-file',
    help="""file .accept for ribotyper.
    It should be one of the followong files:
    `data/archaea.accept`
    `data/bacteria.accept`""",
    required=True
)

# Parameters
parser.add_argument(
    '--ribotyper-threads',
    help='number of threads for ribotyper to use',
    required=False,
    default=1
)

args = parser.parse_args()


# == Import them now ==
import re
import sys
import subprocess as sp
from typing import List

import pandas as pd
from Bio import SeqIO

import src.rg_tools_IO as rgIO


# For convenience
fasta_seqs_fpath = os.path.abspath(args.in_fasta_file)
if not args.prev_short_out_tsv is None \
   and not args.prev_long_out_tsv is None:
    cache_mode = True
    prev_short_out_fpath = os.path.abspath(args.prev_short_out_tsv)
    prev_long_out_fpath  = os.path.abspath(args.prev_long_out_tsv)
else:
    cache_mode = False
    prev_short_out_fpath = None
    prev_long_out_fpath  = None
# end if
outdpath = os.path.abspath(args.outdir)
ribotyper_fpath = os.path.abspath(args.ribotyper)
acccept_fpath = os.path.abspath(args.acccept_file)


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, ribotyper_fpath, acccept_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist')
        sys.exit(1)
    # end if
# enb for

# Check if executables are actually executable
if not os.access(ribotyper_fpath, os.X_OK):
    print(f'Error: file `{ribotyper_fpath}` is not executable')
    sys.exit(1)
# end if

# Create output directory if needed
if not os.path.isdir(outdpath):
    try:
        os.makedirs(outdpath)
    except OSError as err:
        print(f'Error: cannot create directory `{outdpath}`')
        sys.exit(1)
    # end try
# end if

# Check if prev files exist
if cache_mode:
    for f in (prev_short_out_fpath, prev_long_out_fpath):
        if not os.path.exists(f):
            print(f'Error: file `{f}` does not exist')
            sys.exit(1)
        # end if
    # end for
# end if

if args.ribotyper_threads is None:
    ribotyper_threads = 1
else:
    try:
        ribotyper_threads = int(args.ribotyper_threads)
        if ribotyper_threads < 1:
            raise ValueError
        # end if
    except ValueError:
        print('Error: `--ribotyper-threads` is invalid: `{}`'.format(args.ribotyper_threads))
        print('  It must be a positive integer number.')
        sys.exit(1)
    # end try
# end if



print(fasta_seqs_fpath)
if cache_mode:
    print(f'Previous .short.out.tsv file: `{prev_short_out_fpath}`')
    print(f'Previous .long.out.tsv file: `{prev_long_out_fpath}`')
# end if
print(ribotyper_fpath)
print(acccept_fpath)
print()


# Header for reformatted .short.out.tsv files
# and .long.out.tsv
SHORT_OUT_COLNAMES = [
    'target',
    'classification',
    'strnd',
    'pass_fail',
    'unexpected_features',
]

LONG_OUT_COLNAMES = [
    'target',
    'pass_fail',
    'length',
    'fm',
    'fam',
    'domain',
    'model',
    'strnd',
    'ht',
    'tscore',
    'bscore',
    's_per_nt',
    'bevalue',
    'tcov',
    'bcov',
    'bfrom',
    'bto',
    'mfrom',
    'mto',
    'scdiff',
    'scd_per_nt',
    'model',
    'tscore',
    'unexpected_features',
]


def reformat_out_file(raw_short_out_fpath: str,
                       out_colnames: List[str],
                       out_short_out_fpath: str) -> None:
    # Read all lines except of those starting with #
    with open(raw_short_out_fpath, 'rt') as raw_short_out_file:
        lines = list(
            map(
                str.strip,
                filter(
                    lambda x: x[0] != '#',
                    raw_short_out_file.readlines()
                )
            )
        )
    # end with

    # Remove multiple spaces with single space
    # Remove first `idx` column
    squeeze_spaces_pattern = re.compile(r' +')
    for i in range(len(lines)):
        l = re.sub(squeeze_spaces_pattern, '\t', lines[i])
        lines[i] = l.partition('\t')[2] # remove first column
    # end for

    # Write result lines to original file
    with open(out_short_out_fpath, 'wt') as out_short_out_file:
        out_short_out_file.write('\t'.join(out_colnames))
        out_short_out_file.write('\n')
        out_short_out_file.write('\n'.join(lines))
        out_short_out_file.write('\n')
    # end with
# end def


def make_query_file(in_seqs_fpath: str,
                    cache_mode: bool,
                    prev_short_out_fpath=None):
    input_seq_records = rgIO.read_and_filter_fasta(in_seqs_fpath)

    # A temp file for sequences for ribotyper to process
    query_seqs_fpath = os.path.join(
        os.path.dirname(outdpath),
        'query_for_ribotyper.fasta'
    )

    if cache_mode:
        print('Loading previous .short.out.tsv file')

        all_curr_seqIDs = set(map(lambda r: r.id, input_seq_records))

        # Read seqIDs from previous .short.out.tsv file
        prev_short_out_df = pd.read_csv(prev_short_out_fpath, sep='\t')
        prev_seqIDs = set(prev_short_out_df['target'])

        cached_seqIDs = all_curr_seqIDs & prev_seqIDs

        # Find seqIDs of sequences to be processed now by ribotyper
        seqIDs_for_ribotyper = all_curr_seqIDs - cached_seqIDs
        seq_records_for_ribotyper = tuple(
            filter(
                lambda r: r.id in seqIDs_for_ribotyper,
                input_seq_records
            )
        )

        print(
            '{}/{} sequences are cached in the previous .short.out.tsv file' \
                .format(
                    len(input_seq_records) - len(seq_records_for_ribotyper),
                    len(input_seq_records)
                )
        )
        print(
            '{} sequences left to be processed by ribotyper' \
                .format(len(seq_records_for_ribotyper))
        )

        query_seq_records = seq_records_for_ribotyper
        cached_short_out_df = prev_short_out_df.query('target in @cached_seqIDs')
    else:
        query_seq_records = input_seq_records
        cached_short_out_df = None
    # end if

    # Write unprocessed sequences to temporary fasta file
    with open(query_seqs_fpath, 'wt') as tmpfile:
        SeqIO.write(
            query_seq_records,
            tmpfile,
            'fasta'
        )
    # end with

    return query_seqs_fpath, cached_short_out_df
# end def


def run_ribotyper(query_seqs_fpath: str,
                  ribotyper_fpath: str,
                  acccept_fpath: str,
                  outdpath: str,
                  cache_mode: bool,
                  ribotyper_threads: int):
    all_seqs_are_cached = cache_mode \
                          and count_seqs_fasta(query_seqs_fpath) == 0

    if not all_seqs_are_cached:
        # Configure command for ribotyper
        command = ' '.join([
            ribotyper_fpath,
            '-f',
            '-n {}'.format(ribotyper_threads),
            '--minusfail',
            '--scfail',
            '--covfail',
            '--inaccept {}'.format(acccept_fpath),
            query_seqs_fpath,
            outdpath
        ])

        print('Running ribotyper command:')
        print(command)
        print('Ribotyper is running silently. Please, wait')

        pipe = sp.Popen(
            command,
            shell=True,
            stdout=sp.PIPE,
            stderr=sp.PIPE,
            encoding='utf-8'
        )
        out, err = pipe.communicate()

        if pipe.returncode != 0:
            print('Error: ribotyper exited with a non-zero status!')
            print(err)
            sys.exit(1)
        # end if
    else:
        print('All sequence records are already processed.')
        print('The script won\'t run ribotyper')
    # end if

    # Remove temporary fasta file
    try:
        os.unlink(query_seqs_fpath)
    except OSError as err:
        print(f'Cannot remove temporary file `{query_seqs_fpath}`')
        print(err)
        print('Continue anyway\n')
    # end try
# end def


def count_seqs_fasta(fpath):
    return len(
        tuple(
            SeqIO.parse(fpath, 'fasta')
        )
    )
# end def


def get_cached_out_long_df(prev_long_out_fpath, cached_short_out_df):
    cached_seqIDs = set(
        cached_short_out_df['target']
    )
    prev_long_df = pd.read_csv(
        prev_long_out_fpath,
        sep='\t',
        dtypes={
            'target': str,
            'pass_fail': str,
            'length': np.uint16,
            'fm': np.uint16,
            'fam': str,
            'domain': str,
            'model': str,
            'strnd': str,
            'ht': np.uint16,
            'tscore': np.float32,
            'bscore': np.float32,
            's_per_nt': np.float32,
            'bevalue': float,
            'tcov': np.float32,
            'bcov': np.float32,
            'bfrom': np.uint16,
            'bto': np.uint16,
            'mfrom': np.uint16,
            'mto': np.uint16,
            'scdiff': np.float32,
            'scd_per_nt': np.float32,
            'model': str,
            'tscore': np.float32,
            'unexpected_features': str,
        }
    )
    cached_long_df = prev_long_df.query('target in @cached_seqIDs').copy()
    return cached_long_df
# end def


# == Proceed ==

# Configure paths to output files
raw_short_out_fpath = os.path.join(outdpath, 'ribotyper_out.ribotyper.short.out')
final_short_out_fpath = raw_short_out_fpath + '.tsv'
raw_long_out_fpath = os.path.join(outdpath, 'ribotyper_out.ribotyper.long.out')
final_long_out_fpath = raw_long_out_fpath + '.tsv'

query_seqs_fpath, cached_short_out_df = make_query_file(
    fasta_seqs_fpath,
    cache_mode,
    prev_short_out_fpath
)

run_ribotyper(
    query_seqs_fpath,
    ribotyper_fpath,
    acccept_fpath,
    outdpath,
    cache_mode,
    ribotyper_threads
)

# Reformat output .short.out.tsv table
print(
    'Reformatting output file `{}` -> `{}`' \
        .format(raw_short_out_fpath, final_short_out_fpath)
)
reformat_out_file(raw_short_out_fpath, SHORT_OUT_COLNAMES, final_short_out_fpath)
print('Done')
# Reformat output .long.out.tsv table
print(
    'Reformatting output file `{}` -> `{}`' \
        .format(raw_long_out_fpath, final_long_out_fpath)
)
reformat_out_file(raw_long_out_fpath, LONG_OUT_COLNAMES, final_long_out_fpath)
print('Done')

if cache_mode:
    # Add cached .short.out.tsv dataframe rows
    cached_short_out_df.to_csv(
        final_short_out_fpath,
        sep='\t',
        index=False,
        header=False,
        mode='a',
        na_rep='NA'
    )

    # Add cached .long.out.tsv dataframe rows
    cached_long_out_df = get_cached_out_long_df(
        prev_long_out_fpath,
        cached_short_out_df
    )
    cached_long_out_df.to_csv(
        final_long_out_fpath,
        sep='\t',
        index=False,
        header=False,
        mode='a',
        na_rep='NA'
    )
# end if

print('\nCompleted!')
print(final_short_out_fpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
