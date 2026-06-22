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
# 2. `--prev-yet-unasm-asm-accs` -- a file `asm_accs_yet_unasm.txt.gz` from the previous RiboGrove release.
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
    '--prev-yet-unasm-asm-accs',
    help='a file `asm_accs_yet_unasm.txt.gz` from the previous RiboGrove release',
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

import polars as pl
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
if args.prev_asm_sum_final is None or args.prev_yet_unasm_asm_accs is None:
    cache_mode = False
    cache_asm_sum_fpath = None
    cache_yet_unasm_asm_accs_fpath = None
else:
    cache_mode = True
    cache_asm_sum_fpath = os.path.abspath(args.prev_asm_sum_final)
    cache_yet_unasm_asm_accs_fpath = os.path.abspath(args.prev_yet_unasm_asm_accs)
    for fpath in (cache_asm_sum_fpath, cache_yet_unasm_asm_accs_fpath):
        if not os.path.exists(fpath):
            print(f'Error: file `{fpath}` does not exist!')
            sys.exit(1)
        # end if
    # end for
# end if


# Check existance of the input files
fpaths_to_check = (infpath, replicon_map_fpath, release_catalog_fpath)
if cache_mode:
    fpaths_to_check = fpaths_to_check + (cache_asm_sum_fpath, cache_yet_unasm_asm_accs_fpath)
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
    print(cache_yet_unasm_asm_accs_fpath)
# end if
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
        '{} -- Removed assemblies having at least one "unplaced-scaffold"' \
            .format(get_time())
    )
    print('  {:,} genomes are removed'.format(len(asm_accs_to_rm)))

    return filt_asm_sum_df
# end def


def remove_yet_unassembled(asm_sum_df,
                           genomes_dirpath,
                           cache_passing_asm_accs,
                           cache_yet_unasm_asm_accs,
                           yet_unasm_outfpath):
    yet_unasm_asm_accs = set()

    sys.stdout.write('{} -- 0/{} genomes checked'.format(get_time(), asm_sum_df.shape[0]))
    sys.stdout.flush()

    for i, row in enumerate(asm_sum_df.to_dicts(), 1):
        asm_acc = row['asm_acc']

        # Try to hit cache
        if asm_acc in cache_passing_asm_accs:
            continue
        # end if
        if asm_acc in cache_yet_unasm_asm_accs:
            yet_unasm_asm_accs.add(asm_acc)
            continue
        # end if

        seqannot_fpath = get_genome_seqannot_fpath(asm_acc, genomes_dirpath)
        if check_if_genome_is_yet_unasm(seqannot_fpath):
            yet_unasm_asm_accs.add(asm_acc)
        # end if
        sys.stdout.write(
            '\r{} -- {}/{} genomes checked' \
                .format(get_time(), i, asm_sum_df.shape[0])
        )
        sys.stdout.flush()
    # end for

    sys.stdout.write(
        '\r{} -- {}/{} genomes checked' \
            .format(get_time(), i, asm_sum_df.shape[0])
    )
    sys.stdout.flush()

    filt_asm_sum_df = asm_sum_df.filter(
        ~pl.col('asm_acc').is_in(yet_unasm_asm_accs)
    )

    # Save "yet unassembled" accessions
    with gzip.open(yet_unasm_outfpath, 'wt') as yet_unasm_file:
        for asm_acc in yet_unasm_asm_accs:
            yet_unasm_file.write('{}\n'.format(asm_acc))
        # end for
    # end with

    print(
        '\n{} -- Removed 1) genomes with sequences containing NNN; and 2) genomes with "map unlocalized" sequences' \
            .format(get_time())
    )
    print(
        '  {:,} genomes have been removed' \
            .format(len(yet_unasm_asm_accs))
    )
    print('  (their assembly accessions are written to `{}`)'.format(yet_unasm_outfpath))
    return filt_asm_sum_df
# end def

def check_if_genome_is_yet_unasm(seqannot_fpath):
    with gzip.open(seqannot_fpath, 'rt') as gbk_file:
        seq_records = SeqIO.parse(gbk_file, 'gb')
        for seq_record in seq_records:
            if seq_contains_NNN(seq_record):
                return True
            # end if
            if 'MAP UNLOCALIZED' in seq_record.description.upper():
                return True
            # end if
        # end for
    # end with

    return False
# end def

def seq_contains_NNN(seq_record):
    return 'NNN' in seq_record.seq
# end def


def get_cache_passing_asm_accs(cache_asm_sum_fpath):
    cache_asm_sum_df = rgIO.read_ass_sum_file(cache_asm_sum_fpath)
    return set(cache_asm_sum_df['asm_acc'])
# end def

def get_cache_yet_unasm_asm_accs(cache_yet_unasm_asm_accs_fpath):
    with gzip.open(cache_yet_unasm_asm_accs_fpath, 'rt') as infile:
        yet_unasm_asm_accs = frozenset(
            map(
                lambda x: x.strip(),
                infile.readlines()
            )
        )
    # end with
    return yet_unasm_asm_accs
# end def



# == Proceed ==

# Read input
in_asm_sum_df = rgIO.read_ass_sum_file(infpath)

replicon_map_df = pl.read_csv(
    replicon_map_fpath,
    separator='\t'
)

if cache_mode:
    cache_passing_asm_accs = get_cache_passing_asm_accs(cache_asm_sum_fpath)
    cache_yet_unasm_asm_accs = get_cache_yet_unasm_asm_accs(cache_yet_unasm_asm_accs_fpath)
else:
    cache_passing_asm_accs = frozenset()
    cache_yet_unasm_asm_accs    = frozenset()
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
print(
    '  {:,} genomes are retained'.format(filt_asm_sum_df.shape[0])
)
print()


# == 2. Remove assemblies having "unplaced-scaffolds" ==

print('2. Removing assemblies having at least one "unplaced-scaffold"')
print('{} -- Start'.format(get_time()))
filt_asm_sum_df = remove_unplaced_scaffolds(
    filt_asm_sum_df,
    replicon_map_df
)
print(
    '  {:,} genomes are retained'.format(filt_asm_sum_df.shape[0])
)
print()


# == 3. Remove genomes with NNN in their sequences ==

print('3. Removing 1) genomes with sequences containing NNN; 2) genomes with "map unlocalized" sequences')
print('{} -- Start'.format(get_time()))
yet_unasm_outfpath = os.path.join(
    os.path.dirname(outfpath),
    'asm_accs_yet_unasm.txt.gz'
)
filt_asm_sum_df = remove_yet_unassembled(
    filt_asm_sum_df,
    genomes_dirpath,
    cache_passing_asm_accs,
    cache_yet_unasm_asm_accs,
    yet_unasm_outfpath
)
print(
    '  {:,} genomes are retained for further work' \
        .format(filt_asm_sum_df.shape[0])
)


# == Write output ==

filt_asm_sum_df.write_csv(
    outfpath,
    separator='\t',
    include_header=True,
    null_value='NA'
)

print('\n{} -- Completed!'.format(get_time()))
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
