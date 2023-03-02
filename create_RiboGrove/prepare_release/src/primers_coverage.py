
import os
import gzip
from functools import partial

import pandas as pd
from Bio import SeqIO

from src.formatting import format_int_number, format_float_number


_PRIMER_KEYS = (
    '1115F-1492R',
    '27F-1492R',
    '27F-338R',
    '27F-534R',
    '341F-785R',
    '341F-944R',
    '515F-1100R',
    '515F-944R',
    '784F-1100R',
    '784F-1193R',
    '939F-1193R',
    '939F-1378R',
)


def make_ribogrove_primer_coverage_df(input_dirpath,
                                      gene_stats_df):
    path_to_raw_tables = _get_path_to_raw_tables(input_dirpath)
    per_phylum_genome_count_df = _count_genomes_per_phylum(
        gene_stats_df
    )
    primers_coverage_df = _count_coverage_per_phylum(
        path_to_raw_tables,
        per_phylum_genome_count_df,
        gene_stats_df
    )

    print(primers_coverage_df)
    return primers_coverage_df
# end def


def _get_path_to_raw_tables(input_dirpath):
    primer_keys_to_final_names = {
        '1115F-1492R': '1115F-1492R.tsv',
        '27F-1492R'  : '27F-1492R.tsv',
        '27F-338R'   : '27F-338R.tsv',
        '27F-534R'   : '27F-534R.tsv',
        '341F-785R'  : '341F-785R.tsv',
        '341F-944R'  : '341F-944R.tsv',
        '515F-1100R' : '515F-1100R.tsv',
        '515F-944R'  : '515F-944R.tsv',
        '784F-1100R' : '784F-1100R.tsv',
        '784F-1193R' : '784F-1193R.tsv',
        '939F-1193R' : '939F-1193R.tsv',
        '939F-1378R' : '939F-1378R.tsv',
    }

    primer_keys_to_fpaths = {
        k : os.path.join(input_dirpath, v) for k, v in primer_keys_to_final_names.items()
    }

    for f in primer_keys_to_fpaths.values():
        if not os.path.exists(f):
            err_msg = 'Error: file `{}` does not exist'.format(f)
            raise OSError(err_msg)
        # end def
    # end for

    return primer_keys_to_fpaths
# end def


def _count_genomes_per_phylum(gene_stats_df):

    bacterial_seqIDs = set(
        gene_stats_df[
            gene_stats_df['Domain'] == 'Bacteria'
        ]['seqID']
    )

    seqID_df = pd.DataFrame({'seqID': tuple(bacterial_seqIDs)})
    del bacterial_seqIDs

    per_phylum_genome_count_df = seqID_df \
        .merge(
            gene_stats_df[['seqID', 'asm_acc', 'Phylum']],
            on='seqID',
            how='left'
        ) \
        .groupby('Phylum', as_index=False) \
        .agg({'asm_acc': lambda x: x.nunique()}) \
        .rename(columns={'asm_acc': 'num_genomes'}) \
        .sort_values(by='num_genomes', ascending=False)


    return per_phylum_genome_count_df
# end def


def _count_coverage_per_phylum(path_to_raw_tables,
                               per_phylum_genome_count_df,
                               gene_stats_df):
    primers_coverage_df = per_phylum_genome_count_df \
        .copy() \
        .sort_values(by='Phylum')

    for primer_key, raw_table_fpath in path_to_raw_tables.items():

        print(primer_key)

        tmp_df = pd.read_csv(raw_table_fpath, sep='\t')

        tmp_df = tmp_df.merge(
            gene_stats_df[['asm_acc', 'Phylum']],
            on='asm_acc',
            how='outer'
        ) \
        .drop_duplicates() \
        .dropna(subset=['seqID'])

        primers_coverage_df['num_revealed_genomes'] = tmp_df \
            .groupby(['Phylum'], as_index=False) \
            .agg({'asm_acc': lambda x: x.nunique()}) \
            .rename(columns={'asm_acc': 'num_revealed_genomes'}) \
            .merge(
                per_phylum_genome_count_df,
                on='Phylum',
                how='right'
            ).sort_values(by='Phylum') \
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

    return primers_coverage_df
# end def


def format_primer_coverage_df(primer_coverage_df,
                              thousand_separator,
                              decimal_separator):

    curr_format_int_number = partial(
        format_int_number,
        thousand_separator=thousand_separator
    )

    curr_format_float_number = partial(
        format_float_number,
        thousand_separator=thousand_separator,
        decimal_separator=decimal_separator,
        digits=2
    )

    fmt_primer_df = _sort_rows_and_columns(primer_coverage_df)

    fmt_primer_df['Phylum'] = fmt_primer_df['Phylum'] \
        .map(_format_phylum_name)

    fmt_primer_df['num_genomes'] = fmt_primer_df['num_genomes'] \
        .map(curr_format_int_number)

    for primer_key in _PRIMER_KEYS:
        fmt_primer_df[primer_key] = fmt_primer_df[primer_key] \
            .map(curr_format_float_number)
    # end for


    return fmt_primer_df
# end def


def _sort_rows_and_columns(primer_df):
    ordered_primer_names = [
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

    fmt_primer_df = primer_df[
        ['Phylum', 'num_genomes'] + ordered_primer_names
    ].sort_values(by='num_genomes', ascending=False)

    return fmt_primer_df.copy()
# end def


def _format_phylum_name(phylum_name):
    return phylum_name.replace('Candidatus ', 'Ca. ')
# end def
