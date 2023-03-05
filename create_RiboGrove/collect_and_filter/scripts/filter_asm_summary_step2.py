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

### Cache files:
# 1. `--prev-asm-sum-final` -- an assembly summary file after the 2nd step of filtering
#   from the previous RiboGrove release.
#   Optional.
# 2. `--prev-NNN-asm-accs` -- a file `asm_accs_NNN.txt.gz` from the previous RiboGrove release.
#   Optional.

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

parser.add_argument(
    '-g',
    '--genomes-dir',
    help='directory that contains downloaded gbff.gz files',
    required=True
)

# Cache files

parser.add_argument(
    '--prev-asm-sum-final',
    help="""an assembly summary file after the 2nd step of filtering
    from the previous RiboGrove release""",
    required=False
)

parser.add_argument(
    '--prev-NNN-asm-accs',
    help='a file `asm_accs_NNN.txt.gz` from the previous RiboGrove release',
    required=False
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

import pandas as pd
from Bio import SeqIO

import src.rg_tools_IO as rgIO
from src.rg_tools_time import get_time
from src.file_navigation import get_genome_seqannot_fpath


infpath = os.path.realpath(args.in_asm_sum)
replicon_map_fpath = os.path.realpath(args.replicon_map)
release_catalog_fpath = os.path.realpath(args.refseq_catalog)
genomes_dirpath = os.path.realpath(args.genomes_dir)
outfpath = os.path.realpath(args.out_asm_sum)

# Parse "cache" arguments
if args.prev_asm_sum_final is None or args.prev_NNN_asm_accs is None:
    cache_mode = False
    cache_asm_sum_fpath      = None
    cache_NNN_asm_accs_fpath = None
else:
    cache_mode = True
    cache_asm_sum_fpath      = os.path.abspath(args.prev_asm_sum_final)
    cache_NNN_asm_accs_fpath = os.path.abspath(args.prev_NNN_asm_accs)
    for fpath in (cache_asm_sum_fpath, cache_NNN_asm_accs_fpath):
        if not os.path.exists(fpath):
            print(f'Error: file `{fpath}` does not exist!')
            sys.exit(1)
        # end if
    # end for
# end if


# Check existance of the input files
fpaths_to_check = (infpath, replicon_map_fpath, release_catalog_fpath)
if cache_mode:
    fpaths_to_check = fpaths_to_check + (cache_asm_sum_fpath, cache_NNN_asm_accs_fpath)
# end if
for fpath in fpaths_to_check:
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# end for
del fpaths_to_check

if not os.path.isdir(genomes_dirpath):
    print('Error: directory `{}` does not exist'.format(genomes_dirpath))
    sys.exit(1)
# end if

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
if cache_mode:
    print(cache_asm_sum_fpath)
    print(cache_NNN_asm_accs_fpath)
# end if
print()


def remove_nonrelease_genomes(in_asm_sum_df,
                              replicon_map_df,
                              release_catalog_fpath,
                              nonrelease_outfpath):
    all_seq_accs = set(replicon_map_df['seq_acc'])
    print('Reading large release-catalog file silently...')
    curr_release_accs = get_curr_release_seq_accs(release_catalog_fpath)

    print('Filtering...')
    nonrelease_seq_accs = all_seq_accs - curr_release_accs
    nonrelease_asm_accs = set(
        replicon_map_df.query('seq_acc in @nonrelease_seq_accs')['asm_acc']
    )
    # Filter remaining sequences
    filt_asm_sum_df = in_asm_sum_df.query('not asm_acc in @nonrelease_asm_accs')

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
    print(
        '  {:,} genomes are retained for further work'.format(filt_asm_sum_df.shape[0])
    )

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


def remove_NNN_genomes(asm_sum_df,
                       genomes_dirpath,
                       cache_nonNNN_asm_accs,
                       cache_NNN_asm_accs,
                       NNN_outfpath):
    NNN_asm_accs = set()

    sys.stdout.write('{} -- 0/{} genomes checked'.format(get_time(), asm_sum_df.shape[0]))
    sys.stdout.flush()

    for i, row in asm_sum_df.iterrows():
        asm_acc = row['asm_acc']

        # Try to hit cache
        if asm_acc in cache_nonNNN_asm_accs:
            continue
        # end if
        if asm_acc in cache_NNN_asm_accs:
            NNN_asm_accs.add(asm_acc)
            continue
        # end if

        seqannot_fpath = get_genome_seqannot_fpath(asm_acc, genomes_dirpath)
        if check_if_genome_is_NNN(seqannot_fpath):
            NNN_asm_accs.add(asm_acc)
        # end if
        sys.stdout.write(
            '\r{} -- {}/{} genomes checked' \
                .format(get_time(), i+1, asm_sum_df.shape[0])
        )
        sys.stdout.flush()
    # end for

    sys.stdout.write(
        '\r{} -- {}/{} genomes checked' \
            .format(get_time(), i+1, asm_sum_df.shape[0])
    )
    sys.stdout.flush()

    filt_asm_sum_df = asm_sum_df.query('not asm_acc in @NNN_asm_accs')

    # Save "NNN" accessions
    with gzip.open(NNN_outfpath, 'wt') as NNN_file:
        for asm_acc in NNN_asm_accs:
            NNN_file.write('{}\n'.format(asm_acc))
        # end for
    # end with

    print(
        '\n{} -- Removed genomes with sequences containing NNN' \
            .format(get_time())
    )
    print(
        '  {:,} genomes have sequences containing NNN' \
            .format(len(NNN_asm_accs))
    )
    print('  (their assembly accessions are written to `{}`)'.format(NNN_outfpath))
    print(
        '  {:,} genomes are retained for further work' \
            .format(filt_asm_sum_df.shape[0])
    )
    return filt_asm_sum_df
# end def

def check_if_genome_is_NNN(seqannot_fpath):
    with gzip.open(seqannot_fpath, 'rt') as gbk_file:
        seq_records = tuple(SeqIO.parse(gbk_file, 'gb'))
    # end with

    for seq_record in seq_records:
        if seq_contains_NNN(seq_record):
            return True
        # end if
    # end for
    return False
# end def

def seq_contains_NNN(seq_record):
    return 'NNN' in seq_record.seq
# end def


def get_cache_nonNNN_asm_accs(cache_asm_sum_fpath):
    cache_asm_sum_df = rgIO.read_ass_sum_file(cache_asm_sum_fpath)
    return set(cache_asm_sum_df['asm_acc'])
# end def

def get_cache_NNN_asm_accs(cache_NNN_asm_accs_fpath):
    with gzip.open(cache_NNN_asm_accs_fpath, 'rt') as infile:
        NNN_asm_accs = set(
            map(
                lambda x: x.strip(),
                infile.readlines()
            )
        )
    # end with
    return NNN_asm_accs
# end def



# == Proceed ==

# Read input
in_asm_sum_df = rgIO.read_ass_sum_file(infpath)

replicon_map_df = pd.read_csv(
    replicon_map_fpath,
    sep='\t'
)

if cache_mode:
    cache_nonNNN_asm_accs = get_cache_nonNNN_asm_accs(cache_asm_sum_fpath)
    cache_NNN_asm_accs    = get_cache_NNN_asm_accs(cache_NNN_asm_accs_fpath)
else:
    cache_nonNNN_asm_accs = set()
    cache_NNN_asm_accs    = set()
# end def


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
print()

# == 2. Remove genomes with NNN in their sequences ==

print('2. Removing genomes with sequences containing NNN')
print('{} -- Start'.format(get_time()))
NNN_outfpath = os.path.join(
    os.path.dirname(outfpath),
    'asm_accs_NNN.txt.gz'
)
filt_asm_sum_df = remove_NNN_genomes(
    filt_asm_sum_df,
    genomes_dirpath,
    cache_nonNNN_asm_accs,
    cache_NNN_asm_accs,
    NNN_outfpath
)


# == Write output ==

with gzip.open(outfpath, 'wt') as outfile:
    filt_asm_sum_df.to_csv(
        outfile,
        sep='\t',
        index=False,
        header=True,
        na_rep='NA',
        encoding='utf-8',
    )
# end with

print('\n{} -- Completed!'.format(get_time()))
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
