
from functools import partial

import pandas as pd

from src.formatting import format_int_number, format_float_number


def make_ribogrove_copy_number_df(gene_stats_df):

    series_nunique = lambda x: x.nunique()

    by_genome_copy_number_df = gene_stats_df.groupby('ass_id', as_index=False) \
        .agg({'seqID': series_nunique}) \
        .rename(columns={'seqID': 'copy_number'}) \
        .merge(
            gene_stats_df[['ass_id', 'species']].drop_duplicates(),
            on='ass_id',
            how='left'
        )

    by_species_copy_number_df = by_genome_copy_number_df.groupby('species', as_index=False) \
        .agg({'copy_number': 'median'})

    by_species_copy_number_df['copy_number'] = by_species_copy_number_df['copy_number'].map(int)

    ribogrove_copy_number_df = by_species_copy_number_df.groupby('copy_number', as_index=False) \
        .agg({'species': series_nunique}) \
        .rename(columns={'species': 'number_of_species'})

    total_species_count = gene_stats_df['species'].nunique()

    ribogrove_copy_number_df['percent_of_species'] = ribogrove_copy_number_df['number_of_species'] \
                                                     / total_species_count \
                                                     * 100

    ribogrove_copy_number_df = ribogrove_copy_number_df.sort_values(
        by='copy_number',
        ascending=True
    )

    print(ribogrove_copy_number_df)

    return ribogrove_copy_number_df
# end def make_ribogrove_copy_number_df


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
    fmt_ribogrove_copy_number_df['number_of_species'] = fmt_ribogrove_copy_number_df['number_of_species'] \
        .map(curr_format_int_number)
    fmt_ribogrove_copy_number_df['percent_of_species'] = fmt_ribogrove_copy_number_df['percent_of_species'] \
        .map(curr_format_float_number)

    return fmt_ribogrove_copy_number_df
# end def format_copy_number_df
