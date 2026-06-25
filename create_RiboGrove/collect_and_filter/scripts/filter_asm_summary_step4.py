#!/usr/bin/env python3

# The script filters Assembly Summary:
#   it removes genomes with sequences containing at least 3 Ns in row.

## Command line arguments

### Input files:
# 1. `-i / --in-asm-sum` -- an assembly summary file after the 1st step of filtering.
#   Mandatory.
# 2. `-g / --genomes-dir` -- a directory where the downlaoded genomes are located.
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
    help='an assembly summary file after the 3rd step of filtering',
    required=True
)

parser.add_argument(
    '-g',
    '--genomes-dir',
    help='directory that contains downloaded gbff.gz files',
    required=True
)

parser.add_argument(
    '-m',
    '--replicon-map',
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
genomes_dirpath = os.path.realpath(args.genomes_dir)
replicon_map_fpath = os.path.realpath(args.replicon_map)
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
fpaths_to_check = (infpath, replicon_map_fpath,)
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
if cache_mode:
    print(cache_asm_sum_fpath)
    print(cache_yet_unasm_asm_accs_fpath)
# end if
print()

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


def make_genome_size_df(replicon_map_fpath: str) -> pl.DataFrame:
    genome_size_df = (
        pl.scan_csv(replicon_map_fpath, separator='\t')
        .group_by('asm_acc')
        .agg(pl.col('seq_len').sum().alias('genome_size'))
        .collect()
    )
    return genome_size_df
# end def



# == Proceed ==

# Read input
in_asm_sum_df = rgIO.read_ass_sum_file(infpath)

if cache_mode:
    cache_passing_asm_accs = get_cache_passing_asm_accs(cache_asm_sum_fpath)
    cache_yet_unasm_asm_accs = get_cache_yet_unasm_asm_accs(cache_yet_unasm_asm_accs_fpath)
else:
    cache_passing_asm_accs = frozenset()
    cache_yet_unasm_asm_accs    = frozenset()
# end def


print('Removing 1) genomes with sequences containing NNN; 2) genomes with "map unlocalized" sequences')
print('{} -- Start'.format(get_time()))
yet_unasm_outfpath = os.path.join(
    os.path.dirname(outfpath),
    'asm_accs_yet_unasm.txt.gz'
)
filt_asm_sum_df = remove_yet_unassembled(
    in_asm_sum_df,
    genomes_dirpath,
    cache_passing_asm_accs,
    cache_yet_unasm_asm_accs,
    yet_unasm_outfpath
)
print(
    '  {:,} genomes are retained for further work' \
        .format(filt_asm_sum_df.shape[0])
)

genome_size_df = make_genome_size_df(replicon_map_fpath)
filt_asm_sum_df = filt_asm_sum_df.join(genome_size_df, on='asm_acc', how='left')


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
