#!/usr/bin/env python3

import os
import sys
import json
import argparse
from functools import reduce
from collections import OrderedDict

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
    '-t',
    '--taxonomy',
    help='file `bacteria/taxonomy/taxonomy.tsv`',
    required=True
)

parser.add_argument(
    '-d',
    '--target-domain',
    help='either `bacteria` or `archaea`',
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
taxonomy_fpath = os.path.abspath(args.taxonomy)
target_domain = args.target_domain.lower().capitalize()
outfpath = os.path.abspath(args.outfile)


if not os.path.isdir(primers_dirpath):
    print('Error: directory `{}` does not exist'.format(primers_dirpath))
    sys.exit(1)
# end if
if not os.path.isfile(taxonomy_fpath):
    print('Error: file `{}` does not exist'.format(taxonomy_fpath))
    sys.exit(1)
# end if
if not target_domain in ('Bacteria', 'Archaea'):
    print('Error: target domain is invalid: `{}`'.format(target_domain))
    print('Valid options: `Bacteria`, `Archaea`')
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

print('Run parameters:')
print('primers_dirpath: `{}`'.format(primers_dirpath))
print('taxonomy_fpath: `{}`'.format(taxonomy_fpath))
print('target_domain: `{}`'.format(target_domain))
print('outfpath: `{}`'.format(outfpath))


# == Constants ==

RANKS = (
    'Kingdom',
    'Phylum',
    'Class',
    'Order',
    'Family',
    'Genus',
    'Species',
)


# == Functions ==

def parse_primer_pairs():
    primers_pairs_fpath = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'collect_and_filter', 'scripts', 'data', 'primers', 'primer_pairs.json'
    )
    with open(primers_pairs_fpath, 'rt') as infile:
        primer_pairs = json.load(infile)
    # end with
    return primer_pairs
# end def


def make_all_primer_pair_dict(primer_pairs):
    all_primer_pair_key_pairs = reduce(
        lambda list_a, list_b: list_a + list_b,
        primer_pairs.values()
    )
    all_primer_pair_dict = OrderedDict()
    for nameF, nameR, v_region_name in all_primer_pair_key_pairs:
        all_primer_pair_dict['{}-{}'.format(nameF, nameR)] = v_region_name
    # end for
    return all_primer_pair_dict
# end def


def make_total_primer_coverage_df(primers_dirpath,
                                  all_primer_pair_dict,
                                  taxonomy_fpath,
                                  target_domain):

    taxonomy_df = pd.read_csv(taxonomy_fpath, sep='\t')
    taxonomy_df = taxonomy_df[taxonomy_df['Domain'] == target_domain]

    paths_to_raw_tables = _get_paths_to_raw_tables(
        primers_dirpath,
        all_primer_pair_dict
    )

    partial_out_dfs = [None] * len(RANKS)

    for i, rank in enumerate(RANKS):
        print('== {} =='.format(rank))
        partial_out_dfs[i] = _make_per_rank_primer_coverage_df(
            primers_dirpath,
            all_primer_pair_dict,
            taxonomy_df,
            paths_to_raw_tables,
            rank
        )
    # end for

    total_cov_df = pd.concat(partial_out_dfs)

    return _rename_sort_columns_final(total_cov_df, all_primer_pair_dict)
# end def


def _make_per_rank_primer_coverage_df(primers_dirpath,
                                      all_primer_pair_dict,
                                      taxonomy_df,
                                      paths_to_raw_tables,
                                      rank='Phylum'):
    per_rank_genome_count_df = _count_genomes_per_rank(
        taxonomy_df,
        rank
    )
    primers_coverage_df = _count_coverage_per_rank(
        paths_to_raw_tables,
        per_rank_genome_count_df,
        taxonomy_df,
        rank
    )

    return _sort_rows_and_columns_per_rank(
        primers_coverage_df,
        all_primer_pair_dict,
        taxonomy_df,
        rank
    )
# end def


def _get_paths_to_raw_tables(primers_dirpath, all_primer_pair_dict):
    primer_keys_to_final_names = {
        key: '{}.tsv'.format(key) for key in all_primer_pair_dict.keys()
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


def _count_genomes_per_rank(taxonomy_df, rank='Phylum'):

    per_rank_genome_count_df = taxonomy_df \
        .groupby(rank, as_index=False) \
        .agg({'asm_acc': lambda x: x.nunique()}) \
        .rename(columns={'asm_acc': 'num_genomes'}) \
        .sort_values(by='num_genomes', ascending=False)

    return per_rank_genome_count_df
# end def


def _count_coverage_per_rank(paths_to_raw_tables,
                             per_rank_genome_count_df,
                             taxonomy_df,
                             rank='Phylum'):
    primers_coverage_df = per_rank_genome_count_df \
        .copy() \
        .sort_values(by=rank)

    for primer_key, raw_table_fpath in paths_to_raw_tables.items():

        sys.stdout.write('{} '.format(primer_key))
        sys.stdout.flush()

        primer_pair_df = pd.read_csv(raw_table_fpath, sep='\t')
        primer_pair_df['asm_acc'] = np.repeat('', primer_pair_df.shape[0])
        primer_pair_df = primer_pair_df.apply(set_asm_acc, axis=1)

        primer_pair_df = primer_pair_df.merge(
            taxonomy_df[['asm_acc', rank]],
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
            .reset_index() \
            .infer_objects(copy=False) \
            .fillna(0) \
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

def set_asm_acc(row):
    # TODO: use parse_asm_acc function in src/ribogrove_seqID.py
    row['asm_acc'] = row['seqID'].partition(':')[0]
    return row
# end def


def _sort_rows_and_columns_per_rank(per_rank_cov_df,
                                    all_primer_pair_dict,
                                    taxonomy_df,
                                    rank='Phylum'):
    out_df = per_rank_cov_df.copy()
    out_df['rank'] = np.repeat(rank, out_df.shape[0]) # Fill the column with "Phylum" strings
    out_df = out_df.merge(
        taxonomy_df[['Domain', rank]],
        on=rank,
        how='left'
    ).drop_duplicates()
    out_df = out_df.rename(columns={rank: 'taxon'}) # "Phylum" -> "taxon"
    output_col_names = ['Domain', 'rank', 'taxon', 'num_genomes'] \
                       + list(all_primer_pair_dict.keys())
    return out_df[output_col_names] \
        .sort_values(by=['Domain', 'num_genomes'], ascending=False)
# end def


def _rename_sort_columns_final(total_cov_df, all_primer_pair_dict):
    rename_columns_map = {
        primer_pair_key: '{}; {} (%)'.format(primer_pair_key, v_region_name)
        for primer_pair_key, v_region_name in all_primer_pair_dict.items()
    }
    rename_columns_map['rank']        = 'Rank'
    rename_columns_map['taxon']       = 'Taxon'
    rename_columns_map['num_genomes'] = 'Number of genomes'
    total_cov_df = total_cov_df.rename(columns=rename_columns_map)
    return total_cov_df
# end def


# === Proceed ===

print('Start making total primer coverage table')

primer_pairs = parse_primer_pairs()
all_primer_pair_dict = make_all_primer_pair_dict(primer_pairs)

out_df = make_total_primer_coverage_df(
    primers_dirpath,
    all_primer_pair_dict,
    taxonomy_fpath,
    target_domain
)

out_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    encoding='utf-8',
    na_rep='NA'
)

print('Completed successfully! Have fun!')
