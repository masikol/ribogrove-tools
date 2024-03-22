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
    '--final-gene-stats',
    help="""A TSV file of final per-gene statistics.
    This is the output of the script `merge_bases_categories_taxonomy.py`""",
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

import numpy as np
import pandas as pd

from src.rg_tools_time import get_time
from src.primers import make_primer_pair_key
from src.file_navigation import primer_pair_key_2_outfpath


# For convenience
final_stats_fpath = os.path.abspath(args.final_gene_stats)
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
if not os.path.exists(final_stats_fpath):
    print(f'Error: file `{final_stats_fpath}` does not exist!')
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


print(final_stats_fpath)
if primers_mode:
    print(primers_dirpath)
# end if
print()


def make_basic_gcn_df(final_stats_fpath):
    stats_df = pd.read_csv(final_stats_fpath, sep='\t')
    basic_gcn_df = stats_df.groupby('asm_acc', as_index=False) \
        .agg({'seqID': lambda x: x.nunique()}) \
        .rename(columns={'seqID': '16S_rRNA_gcn'})
    return basic_gcn_df
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
    
    for nameF, nameR in primer_pairs:
        primer_pair_key = make_primer_pair_key(nameF, nameR)
        primers_anneal_df_fpath = primer_pair_key_2_outfpath(primers_dirpath, primer_pair_key)
        primers_anneal_df = pd.read_csv(primers_anneal_df_fpath, sep='\t')

        primers_gcn_df = primers_anneal_df.groupby('asm_acc', as_index=False) \
            .agg({'seqID': lambda x: x.nunique()}) \
            .rename(columns={'seqID': '16S_rRNA_gcn'}) \
            .merge(fake_total_df, on='asm_acc', how='right')
        primers_gcn_df['16S_rRNA_gcn'] = primers_gcn_df['16S_rRNA_gcn'].fillna(0).map(int)
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


# == Proceed ==

# Basic GCN df
basic_gcn_df = make_basic_gcn_df(final_stats_fpath)
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
