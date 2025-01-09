#!/usr/bin/env python3

## Description

# The script calculates 16S rRNA Gene Copy Numbers per genome.
# Also, it can calculate "primer-wise" GCNs.
# In this mode, the script also tests if PCR primers anneal to gene sequences.
# Thus, the script will give a genome +1 copy only if the gene sequence
#   can form a PCR product with a primer pair.

### Input files

# 1. `-f / --final-gene-stats` -- A TSV file of final per-gene statistics.
#   This is the output of the script `merge_bases_categories_taxonomy.py`
#   Mandatory.
# 2. `-p / --primers-dir` -- A directory with primer annealing modelling results.
#   This is the output of the script `check_primers_mfeprimer.py`.
#   Optional.

### Ouput files
# 1. `-o / --outdir` -- output directory.


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
    '--final-base-counts',
    help="""A TSV file of final base counts.
    This is the output of the script `count_bases.py`""",
    required=True
)

parser.add_argument(
    '-p',
    '--primers-dir',
    help="""A directory with primer annealing modelling results.
    This is the output of the script `check_primers_mfeprimer.py`.""",
    required=False
)

# Output files

parser.add_argument(
    '-o',
    '--outdir',
    help='output directory',
    required=True
)

args = parser.parse_args()


# == Import them now ==
import sys
import json
from functools import reduce

import numpy as np
import pandas as pd

from src.rg_tools_time import get_time
from src.ribogrove_seqID import parse_asm_acc
from src.primers import make_primer_pair_key
from src.file_navigation import primer_pair_key_2_outfpath


# For convenience
base_counts_fpath = os.path.abspath(args.final_base_counts)
outdir_path = os.path.abspath(args.outdir)
if not args.primers_dir is None:
    print('INFO: the script will calculate primer-wise GCNs as well')
    primers_mode = True
    primers_dirpath = os.path.abspath(args.primers_dir)
else:
    primers_mode = False
    primers_dirpath = None
# end if


# Check existance of all input files and dependencies
if not os.path.exists(base_counts_fpath):
    print(f'Error: file `{base_counts_fpath}` does not exist!')
    sys.exit(1)
# end if

if primers_mode and not os.path.isdir(primers_dirpath):
    print(f'Error: directory `{primers_dirpath}` does not exist!')
    sys.exit(1)
# end if


# Create output directory if needed
if not os.path.isdir(outdir_path):
    try:
        os.makedirs(outdir_path)
    except OSError as err:
        print(f'Error: cannot create directory `{outdir_path}`')
        sys.exit(1)
    # end try
# end if


print(base_counts_fpath)
if primers_mode:
    print(primers_dirpath)
# end if
print()


def make_basic_gcn_df(base_counts_fpath):
    base_count_df = pd.read_csv(base_counts_fpath, sep='\t')
    if 'ass_id' in base_count_df.columns:
        # Consistency with RiboGrove releases before 11.217
        base_count_df = base_count_df.rename(columns={'ass_id': 'asm_acc'})
    elif not 'asm_acc' in base_count_df.columns:
        # Consistency with RiboGrove releases after 21.227
        base_count_df['asm_acc'] = np.repeat('', base_count_df.shape[0])
        base_count_df = base_count_df.apply(set_asm_acc, axis=1)
    # end if
    basic_gcn_df = base_count_df.groupby('asm_acc', as_index=False) \
        .agg({'seqID': lambda x: x.nunique()}) \
        .rename(columns={'seqID': '16S_rRNA_gcn'})
    return basic_gcn_df
# end def

def set_asm_acc(row):
    row['asm_acc'] = parse_asm_acc(row['seqID'])
    return row
# end def

def output_gcn_df(df, outfpath):
    df.to_csv(
        outfpath,
        sep='\t',
        encoding='utf-8',
        index=False,
        header=True,
        na_rep='NA'
    )
# end def


def make_primers_gcn_dfs(primer_pairs,
                         primers_dirpath,
                         basic_gcn_df,
                         outdir_path):

    fake_total_df = pd.DataFrame(
        {
            'asm_acc': list(basic_gcn_df['asm_acc'])
        }
    )
    
    for nameF, nameR, _ in primer_pairs:
        primer_pair_key = make_primer_pair_key(nameF, nameR)
        primers_anneal_df_fpath = primer_pair_key_2_outfpath(primers_dirpath, primer_pair_key)
        primers_anneal_df = pd.read_csv(primers_anneal_df_fpath, sep='\t')
        primers_anneal_df['asm_acc'] = np.repeat('', primers_anneal_df.shape[0])
        primers_anneal_df = primers_anneal_df.apply(set_asm_acc, axis=1)

        primers_gcn_df = primers_anneal_df.groupby('asm_acc', as_index=False) \
            .agg({'seqID': lambda x: x.nunique()}) \
            .rename(columns={'seqID': '16S_rRNA_gcn'}) \
            .merge(fake_total_df, on='asm_acc', how='right')
        primers_gcn_df['16S_rRNA_gcn'] = primers_gcn_df['16S_rRNA_gcn'] \
            .infer_objects(copy=False).fillna(0).map(np.uint8) # we use uint8, for 16S GCN is by no means likely to be >255
        primers_gcn_df = primers_gcn_df[
            ['asm_acc', '16S_rRNA_gcn',]
        ]

        primer_gcn_outfpath = os.path.join(
            outdir_path,
            '{}.tsv'.format(primer_pair_key)
        )
        output_gcn_df(primers_gcn_df, primer_gcn_outfpath)
        print('{}: `{}`'.format(primer_pair_key, primer_gcn_outfpath))
    # end for
# end def


# == Get primer data ==

primers_data_dir = os.path.join(
    os.path.dirname(__file__),
    'data', 'primers'
)

primer_pairs_fpath = os.path.join(
    primers_data_dir,
    'primer_pairs.json'
)
with open(primer_pairs_fpath, 'rt') as infile:
    primer_pairs = json.load(infile)
# end with
primer_pairs = reduce(
    lambda list_a, list_b: list_a + list_b,
    primer_pairs.values()
)


# == Proceed ==

# Basic GCN df
basic_gcn_df = make_basic_gcn_df(base_counts_fpath)
basic_gcn_outfpath = os.path.join(outdir_path, '16S_GCNs.tsv')
output_gcn_df(basic_gcn_df, basic_gcn_outfpath)
print('Basic GCNs: `{}`'.format(basic_gcn_outfpath))


if primers_mode:
    # GCNs for different primer pairs
    make_primers_gcn_dfs(
        primer_pairs,
        primers_dirpath,
        basic_gcn_df,
        outdir_path
    )
# end if

print('\n{} -- Completed!'.format(get_time()))
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
