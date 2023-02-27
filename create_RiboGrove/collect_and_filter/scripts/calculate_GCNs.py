#!/usr/bin/env python3


# TODO: add description

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')

import sys
import json
import argparse

import numpy as np
import pandas as pd

from src.rg_tools_time import get_time
from src.primers import make_primer_pair_key
from src.file_navigation import primer_pair_key_2_outfpath


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--final-gene-stats',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-p',
    '--primers-dir',
    help='fasta file of SSU gene sequences',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outdir',
    help='output directory',
    required=True
)

args = parser.parse_args()


# For convenience
final_stats_fpath = os.path.abspath(args.final_gene_stats)
primers_dirpath = os.path.abspath(args.primers_dir)
outdir_path = os.path.abspath(args.outdir)


# Check existance of all input files and dependencies
if not os.path.exists(final_stats_fpath):
    print(f'Error: file `{final_stats_fpath}` does not exist!')
    sys.exit(1)
# end if

if not os.path.isdir(primers_dirpath):
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
print(primers_dirpath)
print()


def make_basic_gcn_df(final_stats_fpath):
    stats_df = pd.read_csv(final_stats_fpath, sep='\t')
    basic_gcn_df = stats_df.groupby('asm_acc', as_index=False) \
        .agg({'seqID': lambda x: x.nunique()}) \
        .rename(columns={'seqID': 'gcn'})
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
            .rename(columns={'seqID': 'gcn'}) \
            .merge(fake_total_df, on='asm_acc', how='right')
        primers_gcn_df['gcn'] = primers_gcn_df['gcn'].fillna(0).map(int)
        # TODO: remove
        # primers_gcn_df = primers_gcn_df.apply(amend_zero_gcn, axis=1)

        primer_gcn_outfpath = os.path.join(
            outdir_path,
            '{}.tsv'.format(primer_pair_key)
        )
        output_gcn_df(primers_gcn_df, primer_gcn_outfpath)
        print('{}: `{}`'.format(primer_pair_key, primer_gcn_outfpath))
    # end for
# end def

# TODO: remove
# def amend_zero_gcn(row):
#     if row['gcn'] == 0:
#         row['gcn'] = int(4e9)
#     # end if
#     return row
# # end def



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
basic_gcn_outfpath = os.path.join(outdir_path, 'basic_GCNs.tsv')
output_gcn_df(basic_gcn_df, basic_gcn_outfpath)
print('Basic GCNs: `{}`'.format(basic_gcn_outfpath))

make_primers_gcn_dfs(
    primer_pairs,
    primers_dirpath,
    basic_gcn_df,
    outdir_path
)

print('\n{} -- Completed!'.format(get_time()))
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
