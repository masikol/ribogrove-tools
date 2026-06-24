#!/usr/bin/env python3

# The script filters downloaded Assembly Summary again:
#   1) removes genomes which don’t belong to the current RefSeq release
#     using the `.catalog` file;
#   2) removes genomes with sequences containing at least 3 Ns in row.

## Command line arguments

### Input files:
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
    help='an assembly summary file after the 2nd step of filtering',
    required=True
)

parser.add_argument(
    '-m',
    '--replicon-map',
    help="""a replicon map file.
    This is the output of the script `make_replicon_map.py`""",
    required=True
)

parser.add_argument(
    '-a',
    '--refseq-catalog',
    help="""A RefSeq "catalog" file of the current release.
This is the file `RefSeq-releaseXXX.catalog.gz` from here:
https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/.
It is better to filter this file with `filter_refseq_catalog.py` before running current script.
""",
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--out-asm-sum',
    help='an assembly summary file after the 2nd filtering step',
    required=True
)

args = parser.parse_args()


# == Import them now ==
import sys
import gzip

import polars as pl
from Bio import SeqIO

import src.rg_tools_IO as rgIO
from src.rg_tools_time import get_time
from src.file_navigation import get_genome_seqannot_fpath


infpath = os.path.realpath(args.in_asm_sum)
replicon_map_fpath = os.path.realpath(args.replicon_map)
release_catalog_fpath = os.path.realpath(args.refseq_catalog)
outfpath = os.path.realpath(args.out_asm_sum)

# Check existance of the input files
fpaths_to_check = (infpath, replicon_map_fpath, release_catalog_fpath,)
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
print(replicon_map_fpath)
print(release_catalog_fpath)
print()


def remove_nonrelease_genomes(in_asm_sum_df,
                              replicon_map_df,
                              release_catalog_fpath,
                              nonrelease_outfpath):
    all_seq_accs = frozenset(replicon_map_df['seq_acc'])
    print('Reading large release-catalog file silently...')
    curr_release_accs = get_curr_release_seq_accs(release_catalog_fpath)

    print('Filtering...')
    nonrelease_seq_accs = all_seq_accs - curr_release_accs

    nonrelease_asm_accs = frozenset(
        replicon_map_df.filter(
            pl.col('seq_acc').is_in(nonrelease_seq_accs)
        )['asm_acc']
    )
    # Filter remaining sequences

    filt_asm_sum_df = in_asm_sum_df.filter(
        ~pl.col('asm_acc').is_in(nonrelease_asm_accs)
    )
    # Save "nonrelease" accessions
    with gzip.open(nonrelease_outfpath, 'wt') as nonrelease_file:
        for asm_acc in nonrelease_asm_accs:
            nonrelease_file.write('{}\n'.format(asm_acc))
        # end for
    # end with

    print(
        '{} -- Removed genomes not belonging to the current RefSeq release' \
            .format(get_time())
    )
    print(
        '  {:,} genomes do not belong to the current RefSeq release' \
            .format(len(nonrelease_seq_accs))
    )
    print('  (their assembly accessions are written to `{}`)'.format(nonrelease_outfpath))

    return filt_asm_sum_df
# end def

def get_curr_release_seq_accs(release_catalog_fpath):
    # Read the catalog file
    if release_catalog_fpath.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    # end if

    curr_release_accs = set()

    with open_func(release_catalog_fpath, 'rt') as catalog_file:
        acc_column_index = 2
        dir_column_index = 3
        separator = '\t'

        for line in catalog_file:
            line_vals = line.split(separator)
            curr_release_accs.add(line_vals[acc_column_index])
        # end for
    # end with

    return curr_release_accs
# end def


def remove_unplaced_scaffolds(filt_asm_sum_df, replicon_map_df):
    asm_accs_to_rm = frozenset(
        replicon_map_df.filter(pl.col('unplaced_scaffold') == 1)['asm_acc']
    )

    filt_asm_sum_df = filt_asm_sum_df.filter(
        ~pl.col('asm_acc').is_in(asm_accs_to_rm)
    )

    print(
        '{} -- Removed assemblies having at least one "unplaced-scaffold/unlocalized-scaffold"' \
            .format(get_time())
    )
    print('  {:,} genomes are removed'.format(len(asm_accs_to_rm)))

    return filt_asm_sum_df
# end def


# == Proceed ==

# Read input
in_asm_sum_df = rgIO.read_ass_sum_file(infpath)

replicon_map_df = pl.read_csv(
    replicon_map_fpath,
    separator='\t'
)

# == 1. Remove genomes which don't belong to the current RefSeq release ==

print('1. Removing genomes which don\'t belong to the current RefSeq release')
print('{} -- Start'.format(get_time()))
nonrelease_outfpath = os.path.join(
    os.path.dirname(outfpath),
    'asm_accs_nonrelease.txt.gz'
)
filt_asm_sum_df = remove_nonrelease_genomes(
    in_asm_sum_df,
    replicon_map_df,
    release_catalog_fpath,
    nonrelease_outfpath
)
print(
    '  {:,} genomes are retained'.format(filt_asm_sum_df.shape[0])
)
print()


# == 2. Remove assemblies having "unplaced-scaffolds" ==

print('2. Removing assemblies having at least one "unplaced-scaffold/unlocalized-scaffold"')
print('{} -- Start'.format(get_time()))
filt_asm_sum_df = remove_unplaced_scaffolds(
    filt_asm_sum_df,
    replicon_map_df
)
print(
    '  {:,} genomes are retained'.format(filt_asm_sum_df.shape[0])
)
print()


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
