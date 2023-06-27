#!/usr/bin/env python3

import os
import sys
import argparse

import numpy as np
import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input data

parser.add_argument(
    '-p',
    '--primers-dir',
    help='directory `bacteria/primers_coverage`',
    required=True
)

parser.add_argument(
    '-s',
    '--gene-stats',
    help='file `bacteria/gene_stats/per_gene_stats.tsv`',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='output TSV file',
    required=True
)

args = parser.parse_args()

primers_dirpath = os.path.abspath(args.primers_dir)
gene_stats_fpath = os.path.abspath(args.gene_stats)
outfpath = os.path.abspath(args.outfile)


if not os.path.isdir(primers_dirpath):
    print('Error: directory `{}` does not exist'.format(primers_dirpath))
    sys.exit(1)
# end if
if not os.path.isfile(gene_stats_fpath):
    print('Error: file `{}` does not exist'.format(gene_stats_fpath))
    sys.exit(1)
# end if
if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print('Error. Cannot create dir `{}`'.format(os.path.dirname(outfpath)))
        sys.exit(1)
    # end try
# end if


# == Constants ==

RANKS = (
    'Phylum',
    'Class',
    'Order',
    'Family',
    'Genus',
    'Species',
)

PRIMER_PAIR_KEYS = [
    # Full
    '27F-1492R',
    # V1-V2
    '27F-338R',
    # V1-V3
    '27F-534R',
    # V3-V4
    '341F-785R',
    # V3-V5
    '341F-944R',
    # V4
    '515F-806R',
    # V4-V5
    '515F-944R',
    # V4-V6
    '515F-1100R',
    # V5-V6
    '784F-1100R',
    # V5-V7
    '784F-1193R',
    # V6-V7
    '939F-1193R',
    # V6-V8
    '939F-1378R',
    # V7-V9
    '1115F-1492R',
]

# In corresponding order to PRIMER_PAIR_KEYS
V_REGION_NAMES = [
    'Full gene',
    'V1-V2',
    'V1-V3',
    'V3-V4',
    'V3-V5',
    'V4',
    'V4-V5',
    'V4-V6',
    'V5-V6',
    'V5-V7',
    'V6-V7',
    'V6-V8',
    'V7-V9',
]

OUT_COL_NAMES = [
    'rank',
    'taxon',
    'num_genomes',
] + PRIMER_PAIR_KEYS


# == Classes ==

class GeneStatsWrapper:

    def __init__(self, gene_stats_fpath):
        self.gene_stats_df = pd.read_csv(gene_stats_fpath, sep='\t') \
            .query('Domain == "Bacteria"')
        bacterial_seqIDs = set(
            self.gene_stats_df['seqID']
        )
        self.seqID_df = pd.DataFrame({'seqID': tuple(bacterial_seqIDs)})
    # end def
# end class



# == Functions ==


def make_total_primer_coverage_df(primers_dirpath, gene_stats_fpath):

    gene_stats_wrapper = GeneStatsWrapper(gene_stats_fpath)

    paths_to_raw_tables = _get_paths_to_raw_tables(primers_dirpath)

    partial_out_dfs = [None] * len(RANKS)

    for i, rank in enumerate(RANKS):
        print('== {} =='.format(rank))
        partial_out_dfs[i] = make_per_rank_primer_coverage_df(
            primers_dirpath,
            gene_stats_wrapper,
            paths_to_raw_tables,
            rank
        )
    # end for

    total_cov_df = pd.concat(partial_out_dfs)

    return _rename_columns_final(total_cov_df)
# end def


def make_per_rank_primer_coverage_df(primers_dirpath,
                                     gene_stats_wrapper,
                                     paths_to_raw_tables,
                                     rank='Phylum'):
    per_rank_genome_count_df = _count_genomes_per_rank(
        gene_stats_wrapper,
        rank
    )
    primers_coverage_df = _count_coverage_per_rank(
        paths_to_raw_tables,
        per_rank_genome_count_df,
        gene_stats_wrapper.gene_stats_df,
        rank
    )
    return _sort_rows_and_columns_per_rank(primers_coverage_df, rank)
# end def


def _get_paths_to_raw_tables(primers_dirpath):
    primer_keys_to_final_names = {
        key: '{}.tsv'.format(key) for key in PRIMER_PAIR_KEYS
    }

    primer_keys_to_fpaths = {
        k : os.path.join(primers_dirpath, v) for k, v in primer_keys_to_final_names.items()
    }

    for f in primer_keys_to_fpaths.values():
        if not os.path.exists(f):
            err_msg = 'Error: file `{}` does not exist'.format(f)
            raise OSError(err_msg)
        # end def
    # end for

    return primer_keys_to_fpaths
# end def


def _count_genomes_per_rank(gene_stats_wrapper, rank='Phylum'):

    per_rank_genome_count_df = gene_stats_wrapper.seqID_df \
        .merge(
            gene_stats_wrapper.gene_stats_df[['seqID', 'asm_acc', rank]],
            on='seqID',
            how='left'
        ) \
        .groupby(rank, as_index=False) \
        .agg({'asm_acc': lambda x: x.nunique()}) \
        .rename(columns={'asm_acc': 'num_genomes'}) \
        .sort_values(by='num_genomes', ascending=False)

    return per_rank_genome_count_df
# end def


def _count_coverage_per_rank(paths_to_raw_tables,
                             per_rank_genome_count_df,
                             gene_stats_df,
                             rank='Phylum'):
    primers_coverage_df = per_rank_genome_count_df \
        .copy() \
        .sort_values(by=rank)

    for primer_key, raw_table_fpath in paths_to_raw_tables.items():

        sys.stdout.write('{} '.format(primer_key))
        sys.stdout.flush()

        primer_pair_df = pd.read_csv(raw_table_fpath, sep='\t')

        primer_pair_df = primer_pair_df.merge(
            gene_stats_df[['asm_acc', rank]],
            on='asm_acc',
            how='outer'
        ) \
        .drop_duplicates() \
        .dropna(subset=['seqID'])

        primers_coverage_df['num_revealed_genomes'] = primer_pair_df \
            .groupby(rank, as_index=False) \
            .agg({'asm_acc': lambda x: x.nunique()}) \
            .rename(columns={'asm_acc': 'num_revealed_genomes'}) \
            .merge(
                per_rank_genome_count_df,
                on=rank,
                how='right'
            ).sort_values(by=rank) \
            .fillna(0) \
            .reset_index() \
            ['num_revealed_genomes']

        primers_coverage_df[primer_key] = primers_coverage_df['num_revealed_genomes'] \
                                          / primers_coverage_df['num_genomes'] \
                                          * 100
        primers_coverage_df = primers_coverage_df.drop(
            ['num_revealed_genomes'],
            axis=1
        )
    # end for

    print()

    return primers_coverage_df
# end def


def _sort_rows_and_columns_per_rank(per_rank_cov_df, rank='Phylum'):
    out_df = per_rank_cov_df.copy()
    out_df['rank'] = np.repeat(rank, out_df.shape[0]) # Fill the column with "Phylum" strings
    out_df = out_df.rename(columns={rank: 'taxon'}) # "Phylum" -> "taxon"
    return out_df[OUT_COL_NAMES] \
        .sort_values(by='num_genomes', ascending=False)
# end def


def _rename_columns_final(total_cov_df):
    rename_columns_map = {
        primer_pair_key: '{}; {} (%)'.format(primer_pair_key, v_region_name)
        for primer_pair_key, v_region_name in zip(PRIMER_PAIR_KEYS, V_REGION_NAMES)
    }
    rename_columns_map['rank']        = 'Rank'
    rename_columns_map['taxon']       = 'Taxon'
    rename_columns_map['num_genomes'] = 'Number of genomes'
    return total_cov_df.rename(columns=rename_columns_map)
# end def


# === Proceed ===

print('Start making total primer coverage table')

out_df = make_total_primer_coverage_df(primers_dirpath, gene_stats_fpath)
out_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    encoding='utf-8',
    na_rep='NA'
)

print('Completed successfully! Have fun!')
