
from functools import partial

import pandas as pd

from src.formatting import format_int_number, format_float_number


def make_ribogrove_copy_number_df(gene_stats_df):
    bacteria_gcn_df = _make_per_species_median_gcn_df(gene_stats_df, 'Bacteria')
    archaea_gcn_df  = _make_per_species_median_gcn_df(gene_stats_df,  'Archaea')

    merged_gcn_df = bacteria_gcn_df.merge(
        archaea_gcn_df,
        on='copy_number',
        how='left'
    ).fillna(0)

    merged_gcn_df.columns = [
        'copy_number',
        'number_of_species_bacteria',
        'percent_of_species_bacteria',
        'number_of_species_archaea',
        'percent_of_species_archaea',
    ]

    int_colnames = [
        'copy_number',
        'number_of_species_bacteria',
        'number_of_species_archaea',
    ]
    for colname in int_colnames:
        merged_gcn_df[colname] = merged_gcn_df[colname].map(int)
    # end for

    print(merged_gcn_df.shape)
    print(merged_gcn_df)

    return merged_gcn_df
# end def make_ribogrove_copy_number_df


def _make_per_species_median_gcn_df(gene_stats_df, domain_name):
    series_nunique = lambda x: x.nunique()

    domain_gene_stats_df = gene_stats_df[gene_stats_df['domain'] == domain_name]

    by_genome_copy_number_df = domain_gene_stats_df.groupby('ass_id', as_index=False) \
        .agg({'seqID': series_nunique}) \
        .rename(columns={'seqID': 'copy_number'}) \
        .merge(
            domain_gene_stats_df[['ass_id', 'species']].drop_duplicates(),
            on='ass_id',
            how='left'
        )

    by_species_copy_number_df = by_genome_copy_number_df.groupby('species', as_index=False) \
        .agg({'copy_number': 'median'})

    by_species_copy_number_df['copy_number'] = by_species_copy_number_df['copy_number'].map(int)

    ribogrove_copy_number_df = by_species_copy_number_df.groupby('copy_number', as_index=False) \
        .agg({'species': series_nunique}) \
        .rename(columns={'species': 'number_of_species'})

    total_species_count = domain_gene_stats_df['species'].nunique()

    ribogrove_copy_number_df['percent_of_species'] = ribogrove_copy_number_df['number_of_species'] \
                                                     / total_species_count \
                                                     * 100

    ribogrove_copy_number_df = ribogrove_copy_number_df.sort_values(
        by='copy_number',
        ascending=True
    )

    return ribogrove_copy_number_df
# end def


def format_copy_number_df(ribogrove_copy_number_df, thousand_separator, decimal_separator):

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

    fmt_ribogrove_copy_number_df = ribogrove_copy_number_df.copy()
    fmt_ribogrove_copy_number_df['copy_number'] = fmt_ribogrove_copy_number_df['copy_number'] \
        .map(curr_format_int_number)
    fmt_ribogrove_copy_number_df['number_of_species_bacteria'] = \
        fmt_ribogrove_copy_number_df['number_of_species_bacteria'].map(curr_format_int_number)
    fmt_ribogrove_copy_number_df['percent_of_species_bacteria'] = \
        fmt_ribogrove_copy_number_df['percent_of_species_bacteria'].map(curr_format_float_number)
    fmt_ribogrove_copy_number_df['number_of_species_archaea'] = \
        fmt_ribogrove_copy_number_df['number_of_species_archaea'].map(curr_format_int_number)
    fmt_ribogrove_copy_number_df['percent_of_species_archaea'] = \
        fmt_ribogrove_copy_number_df['percent_of_species_archaea'].map(curr_format_float_number)

    return fmt_ribogrove_copy_number_df
# end def format_copy_number_df
