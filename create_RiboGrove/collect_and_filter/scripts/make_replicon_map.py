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

import pandas as pd

import src.rg_tools_IO as rgIO
from src.file_navigation import get_asm_report_fpath


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


def make_replicon_map(asm_sum_fpath, genomes_dirpath, prev_repl_map_fpath):
    asm_sum_df = rgIO.read_ass_sum_file(asm_sum_fpath)
    all_accs = set(asm_sum_df['asm_acc'])

    if cache_mode:
        cached_df = load_prev_repl_map(all_accs, prev_repl_map_fpath)
        cached_asm_accs = set(cached_df['asm_acc'])
    # end if

    final_asm_accs, final_seq_accs = list(), list()
    for _, asm_acc in enumerate(all_accs):
        if cache_mode and asm_acc in cached_asm_accs:
            seq_accessions = tuple(
                cached_df[cached_df['asm_acc'] == asm_acc]['seq_acc']
            )
        else:
            seq_accessions = extract_sec_accessions(asm_acc, genomes_dirpath)
        # end if
        for seq_acc in seq_accessions:
            final_asm_accs.append(asm_acc)
            final_seq_accs.append(seq_acc)
        # end for
    # end for

    replicon_map_df = pd.DataFrame(
        {
            'asm_acc': final_asm_accs,
            'seq_acc': final_seq_accs,
        }
    )
    return replicon_map_df
# end def


def load_prev_repl_map(all_curr_accs, prev_repl_map_fpath):
    with gzip.open(prev_repl_map_fpath, 'rt') as infile:
        prev_df = pd.read_csv(infile, sep='\t')
    # end with
    prev_df = prev_df.query('asm_acc in @all_curr_accs')
    return prev_df
# end def

def extract_sec_accessions(asm_acc, genomes_dirpath):
    asm_report_fpath = get_asm_report_fpath(asm_acc, genomes_dirpath)
    with open(asm_report_fpath, 'rt') as infile:
        lines = remove_comment_lines(infile.readlines())
        accessions = tuple(map(get_refseq_accession, lines))
    # end with
    return accessions
# end def

def remove_comment_lines(lines):
    return tuple(
        filter(
            doesnt_start_with_num_sign,
            lines
        )
    )
# end def

def doesnt_start_with_num_sign(line):
    return line[0] != '#'
# end def

def get_refseq_accession(line):
    return line.split('\t')[6]
# end def


def write_output(replicon_map_df, outfpath):
    with gzip.open(outfpath, 'wt') as outfile:
        replicon_map_df.to_csv(
            outfile,
            sep='\t',
            encoding='utf-8',
            index=False,
            header=True,
            na_rep='NA'
        )
    # end with
# end def



# == Proceed ==

replicon_map_df = make_replicon_map(
    asm_sum_fpath,
    genomes_dirpath,
    prev_repl_map_fpath
)
write_output(replicon_map_df, outfpath)


print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
