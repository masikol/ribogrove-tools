#!/usr/bin/env python3

# The script makes a helper file in which Assembly acceession numbers
#   are mapped to corresponding RefSeq accession numbers.
# For example `GCF_000005825.2` is mapped to `NC_013791.2`, `NC_013792.1` and `NC_013793.1`.

## Command line arguments

### Input files:
# 1. `-i / --asm-sum` -- an assembly summary file after the 1st step of filtering.
#   Mandatory.
# 2. `-g / --genomes-dir` -- a directory where the downlaoded genomes are located.
#   It is the output of the script `download_genomes.py`.
#   Mandatory.

### "Cached" files:
# 1. `--prev-replicon-map` -- a replicon map
#   of the previous RiboGrove release (replicon_map.tsv.gz).

### Output files:
# 1. `-o / --out` -- a file in which Assembly acceession numbers
#   are mapped to corresponding RefSeq accession numbers.
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

# Input
parser.add_argument(
    '-i',
    '--asm-sum',
    help='an assembly summary file after the 1st step of filtering',
    required=True
)

parser.add_argument(
    '-g',
    '--genomes-dir',
    help="""A directory where the downlaoded genomes are located.
    It is the output of the script `download_genomes.py`.""",
    required=True
)

# "Chache" files
parser.add_argument(
    '--prev-replicon-map',
    help="""a replicon map of the previous RiboGrove release""",
    required=False
)

# Output
parser.add_argument(
    '-o',
    '--out',
    help="""A file in which Assembly acceession numbers
    are mapped to corresponding RefSeq accession numbers""",
    required=True
)


args = parser.parse_args()


# == Import them now ==
import sys
import gzip
from io import StringIO

import polars as pl
from Bio import SeqIO

import src.rg_tools_IO as rgIO
from src.file_navigation import get_asm_report_fpath, get_genome_seqannot_fpath

UNPLACED_LABELS = frozenset([
    'unplaced-scaffold',
    'unlocalized-scaffold',
])


asm_sum_fpath = os.path.realpath(args.asm_sum)
outfpath = os.path.realpath(args.out)
genomes_dirpath = os.path.realpath(args.genomes_dir)


# Check existance of input file -c/--gi-2-acc-fpath
if not os.path.exists(asm_sum_fpath):
    print(f'Error: file `{asm_sum_fpath}` does not exist')
    sys.exit(1)
# end if

if not os.path.isdir(genomes_dirpath):
    print(f'Error: directory `{genomes_dirpath}` does not exist')
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

cache_mode = not args.prev_replicon_map is None
if cache_mode:
    prev_repl_map_fpath = os.path.abspath(args.prev_replicon_map)
    if not os.path.isfile(prev_repl_map_fpath):
        print('Error!')
        print('File `{}` does not exist'.format(prev_repl_map_fpath))
        sys.exit(1)
    # end if
else:
    prev_repl_map_fpath = None
# end if

print(asm_sum_fpath)
print(genomes_dirpath)
if cache_mode:
    print(prev_repl_map_fpath)
# end if
print()


def make_replicon_map(asm_sum_fpath,
                      genomes_dirpath,
                      prev_repl_map_fpath,
                      outfpath):
    asm_sum_df = rgIO.read_ass_sum_file(asm_sum_fpath)
    all_accs = set(asm_sum_df['asm_acc'])

    if cache_mode:
        cached_df = load_prev_repl_map(all_accs, prev_repl_map_fpath)
        cached_asm_accs = set(cached_df['asm_acc'])
    # end if

    sep = '\t'

    with gzip.open(outfpath, 'wt') as out_handle:
        out_handle.write('{}\n'.format(sep.join([
            'asm_acc',
            'seq_acc',
            'seq_len',
            'unplaced_scaffold',
        ])))

        for asm_acc in all_accs:
            if cache_mode and asm_acc in cached_asm_accs:
                curr_cache_df = cached_df.filter(pl.col('asm_acc') == asm_acc)
                curr_cache_df.write_csv(
                    curr_cache_df,
                    separator='\t',
                    include_header=False,
                    null_value='NA'
                )
                continue
            # end if
            replicon_info = extract_replicon_info(asm_acc, genomes_dirpath)
            for seq_acc, info_dict in replicon_info.items():
                out_handle.write('{}\n'.format(sep.join([
                    asm_acc,
                    seq_acc,
                    str(info_dict['seq_len']),
                    '1' if info_dict['unplaced_scaffold'] == True else '0',
                ])))
            # end for
        # end for
    # end with
# end def


def load_prev_repl_map(all_curr_accs, prev_repl_map_fpath):
    prev_df = (
        pl.scan_csv(prev_repl_map_fpath, separator='\t')
        .filter(p.col('asm_acc').is_in(all_curr_accs))
        .collect()
    )
    prev_df = prev_df
    return prev_df
# end def

def extract_replicon_info(asm_acc, genomes_dirpath):
    asm_report_fpath = get_asm_report_fpath(asm_acc, genomes_dirpath)

    report_df = pl.read_csv(
        asm_report_fpath,
        separator='\t',
        has_header=False,
        comment_prefix='#',
        columns=[1, 6, 8],
        new_columns=[
            'seq_role',
            'refseq_acc',
            'seq_len'
        ]
    )

    seq_lengths = infer_seq_lengths(report_df, asm_acc, genomes_dirpath)
    unplaced_scaffold_dict = {
        row['refseq_acc']: row['seq_role'] in UNPLACED_LABELS for row in report_df.to_dicts()
    }

    result_dict = {
        seq_acc: {
            'seq_len': seq_lengths[seq_acc],
            'unplaced_scaffold': unplaced_scaffold_dict[seq_acc],
        } for seq_acc in seq_lengths.keys()
    }

    return result_dict
# end def

def infer_seq_lengths(report_df, asm_acc, genomes_dirpath):
    all_lengths_are_ints = False
    try:
        _ = tuple(map(int, report_df['seq_len']))
    except:
        pass
    else:
        all_lengths_are_ints = True
    # end try

    some_seq_accs_are_na = any(map(
        lambda seq_acc: seq_acc.upper() == 'NA',
        report_df['refseq_acc']
    ))

    if all_lengths_are_ints and not some_seq_accs_are_na:
        seq_lengths = {
            row['refseq_acc']: int(row['seq_len']) for row in report_df.to_dicts()
        }
    else:
        print(
            'Inferring seq lengths for {} from gbk...'.format(asm_acc),
            end=' '
        )
        seq_lengths = infer_seq_lengths_from_gbk(asm_acc, genomes_dirpath)
        print('  {} bp'.format(';'.join(
            tuple(map(str, seq_lengths.values()))
        )))
    # end if
    return seq_lengths
# end def

def infer_seq_lengths_from_gbk(asm_acc, genomes_dirpath):
    infpath = get_genome_seqannot_fpath(asm_acc, genomes_dirpath)
    with gzip.open(infpath, 'rt') as in_handle:
        seq_records = tuple(SeqIO.parse(in_handle, 'genbank'))
    # end with
    return {
        sr.id: len(sr.seq) for sr in seq_records
    }
# end def


# == Proceed ==

replicon_map_df = make_replicon_map(
    asm_sum_fpath,
    genomes_dirpath,
    prev_repl_map_fpath,
    outfpath
)


print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
