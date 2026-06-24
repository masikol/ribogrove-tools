#!/usr/bin/env python3

# The script filters downloaded Assembly Summary again:
#   1) removes genomes which don’t belong to the current RefSeq release
#     using the `.catalog` file;
#   2) removes genomes with sequences containing at least 3 Ns in row.

## Command line arguments

### Input files:
# TODO: update
# 1. `-i / --in-asm-sum` -- an assembly summary file after the 1st step of filtering.
#   Mandatory.
# 2. `-m / --replicon-map` -- a replicon map file.
#   This is the output of the script `make_replicon_map.py`
#   Mandatory.
# 3. `-a / --refseq-catalog` -- A RefSeq "catalog" file of the current release.
#   This is the file `RefSeq-releaseXXX.catalog.gz` from here:
#   https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/.
#   It is better to filter this file with `filter_refseq_catalog.py` before running current script.
#   Mandatory.
# 4. `-g / --genomes-dir` -- a directory where the downlaoded genomes are located.
#   It is the output of the script `download_genomes.py`.
#   Mandatory.

### Output files:

# 1. `--outfile` -- an assembly summary file after the 2nd filtering step.
#   Mandatory.


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
    '--in-asm-sum',
    help='an assembly summary file after the 1st step of filtering',
    required=True
)

parser.add_argument(
    '-s',
    '--in-sample-types',
    help='file `sample_types.tsv` produced by the script `make_sample_type_table.py`',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--out-asm-sum',
    help='output assembly summary',
    required=True
)

args = parser.parse_args()


# == Import them now ==
import sys
import gzip

import polars as pl

import src.rg_tools_IO as rgIO


infpath = os.path.realpath(args.in_asm_sum)
sample_type_df_fpath = os.path.realpath(args.in_sample_types)
outfpath = os.path.realpath(args.out_asm_sum)

# Check existance of the input files
fpaths_to_check = (infpath, sample_type_df_fpath,)
for fpath in fpaths_to_check:
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# end for
del fpaths_to_check

if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if


print(infpath)
print(sample_type_df_fpath)
print()

SAMPLE_TYPES_TO_REMOVE = frozenset(map(
    str.lower,
    [
        'Metagenome assemble',
        'Metagenome assembly',
        'metagenomic',
        'metagenomic assembiy',
        'metagenomic assembly',
        'Metagenomic assembly',
        'Metagenomic Assembly',
    ]
))


# == Proceed ==

# Read input
in_asm_sum_df = rgIO.read_ass_sum_file(infpath)

sample_type_df = pl.read_csv(
    sample_type_df_fpath,
    separator='\t',
    null_values=['', 'na', 'NA']
).with_columns(
    pl.col('sample_type').str.to_lowercase()
)

in_asm_sum_df = in_asm_sum_df.join(
    sample_type_df,
    on='asm_acc',
    how='left'
)

num_genomes_before = in_asm_sum_df['asm_acc'].n_unique()
print('Number of input genomes: {:,}'.format(num_genomes_before))
print('Removing genomes that originate from metagenomic samples')
print('{} -- Start'.format(get_time()))

filt_asm_sum_df = in_asm_sum_df.filter(
    ~pl.col('sample_type').is_in(SAMPLE_TYPES_TO_REMOVE, nulls_equal=True)
)

print('{} -- done'.format(get_time()))
num_genomes_after = filt_asm_sum_df['asm_acc'].n_unique()
print(
    '  {:,} genomes are removed'.format(
        num_genomes_before - num_genomes_after
    )
)
print(
    '  {:,} genomes are retained'.format(num_genomes_after)
)
print()

seem_like_metagenomic_df = filt_asm_sum_df.select(
    pl.col('asm_acc', 'sample_type')
).filter(
    pl.col('sample_type').str.contains('metagenom')
)


if seem_like_metagenomic_df.height > 0:
    print('\nWARNING!')
    print('{} genomes still seem to originate from metagenomic samples,'.format(
        seem_like_metagenomic_df.height
    ))
    print('  even though their `sample_type` do not appear in `SAMPLE_TYPES_TO_REMOVE`')
    print('Here they are:')
    for i, row in enumerate(seem_like_metagenomic_df.to_dicts(), 1):
        print('  {}. {} -- `{}`'.format(i, row['asm_acc'], row['sample_type']))
    # end for
    print('And here are `SAMPLE_TYPES_TO_REMOVE`:')
    print('  {}'.format(str(SAMPLE_TYPES_TO_REMOVE)))
    sys.exit(1)
# end if


# == Write output ==

filt_asm_sum_df.write_csv(
    outfpath,
    separator='\t',
    include_header=True,
    null_value='NA',
    compression='gzip'
)

print('\n{} -- Completed!'.format(get_time()))
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
