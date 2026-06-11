
from functools import partial

import polars as pl

from src.util import remove_invalid_species
from src.formatting import format_int_number, format_float_number


def make_ribogrove_copy_number_df(gene_stats_df):

    bacteria_gcn_df = _make_per_species_median_gcn_df(gene_stats_df, 'Bacteria')
    archaea_gcn_df  = _make_per_species_median_gcn_df(gene_stats_df,  'Archaea')

    # TODO: remove pd
    # merged_gcn_df = bacteria_gcn_df.merge(
    #     archaea_gcn_df,
    #     on='copy_number',
    #     how='left'
    # ).fillna(0)
    merged_gcn_df = bacteria_gcn_df.join(
        archaea_gcn_df,
        on='copy_number',
        how='left'
    ).with_columns(pl.col('copy_number').fill_null(0))

    merged_gcn_df.columns = [
        'copy_number',
        'number_of_species_bacteria',
        'percent_of_species_bacteria',
        'number_of_species_archaea',
        'percent_of_species_archaea',
    ]

    # TODO: remove pd
    # int_colnames = [
    #     'copy_number',
    #     'number_of_species_bacteria',
    #     'number_of_species_archaea',
    # ]
    # for colname in int_colnames:
    #     merged_gcn_df[colname] = merged_gcn_df[colname].map(int)
    # # end for
    merged_gcn_df = merged_gcn_df.with_columns(
        pl.col('copy_number').cast(pl.Uint32),
        pl.col('number_of_species_bacteria').cast(pl.Uint32),
        pl.col('number_of_species_archaea').cast(pl.Uint32)
    )

    print(merged_gcn_df.shape)
    print(merged_gcn_df)

    return merged_gcn_df
# end def


def _make_per_species_median_gcn_df(gene_stats_df, domain_name):

    # TODO: remove pd
    # series_nunique = lambda x: x.nunique()

    # TODO: remove pd
    # domain_gene_stats_df = gene_stats_df[gene_stats_df['Domain'] == domain_name] \
    #     .reset_index()
    # tmp_subset_df = domain_gene_stats_df[
    #     is_validlike_species_name(domain_gene_stats_df['Species'])
    # ]
    domain_gene_stats_df = gene_stats_df.filter(pl.col('Domain') == domain_name)
    tmp_subset_df = remove_invalid_species(domain_gene_stats_df)

    # TODO: remove pd
    # by_genome_copy_number_df = tmp_subset_df.groupby('asm_acc', as_index=False) \
    #     .agg({'seqID': series_nunique}) \
    #     .rename(columns={'seqID': 'copy_number'}) \
    #     .merge(
    #         tmp_subset_df[['asm_acc', 'Species']].drop_duplicates(),
    #         on='asm_acc',
    #         how='left'
    #     )
    by_genome_copy_number_df = tmp_subset_df.group_by('asm_acc') \
        .agg(pl.col('seqID').n_unique()) \
        .rename({'seqID': 'copy_number'}) \
        .join(
            tmp_subset_df.select(pl.col('asm_acc', 'Species')).unique(),
            on='asm_acc',
            how='left'
        )

    # TODO: remove pd
    # by_species_copy_number_df = by_genome_copy_number_df.groupby('Species', as_index=False) \
    #     .agg({'copy_number': 'median'})
    by_species_copy_number_df = by_genome_copy_number_df.group_by('Species') \
        .agg(pl.col('copy_number').median())

    # TODO: remove pd
    # by_species_copy_number_df['copy_number'] = by_species_copy_number_df['copy_number'].map(int)
    by_species_copy_number_df = by_species_copy_number_df.with_columns(
        pl.col('copy_number').cast(pl.UInt32)
    )

    # TODO: remove pd
    # ribogrove_copy_number_df = by_species_copy_number_df.groupby('copy_number', as_index=False) \
    #     .agg({'Species': series_nunique}) \
    #     .rename(columns={'Species': 'number_of_species'})
    ribogrove_copy_number_df = by_species_copy_number_df.group_by('copy_number') \
        .agg(pl.col('Species').n_unique()) \
        .rename({'Species': 'number_of_species'})

    # TODO: remove pd
    # total_species_count = tmp_subset_df['Species'].nunique()
    total_species_count = tmp_subset_df['Species'].n_unique()

    # TODO: remove pd
    # ribogrove_copy_number_df['percent_of_species'] = ribogrove_copy_number_df['number_of_species'] \
    #                                                  / total_species_count \
    #                                                  * 100
    ribogrove_copy_number_df = ribogrove_copy_number_df.with_columns(
        (
            pl.col(number_of_species) / total_species_count * 100.0
        ).alias('percent_of_species')
    )

    # TODO: remove pd
    # ribogrove_copy_number_df = ribogrove_copy_number_df.sort_values(
    #     by='copy_number',
    #     ascending=True
    # )
    ribogrove_copy_number_df = ribogrove_copy_number_df.sort(
        by='copy_number',
        descending=False
    )

    return ribogrove_copy_number_df
# end def


def format_copy_number_df(copy_number_df, thousand_separator, decimal_separator):

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

    fmt_copy_number_df = copy_number_df.with_columns(
        pl.col('copy_number').map_elements(curr_format_int_number, return_dtype=pl.String),
        pl.col('number_of_species_bacteria').map_elements(curr_format_int_number, return_dtype=pl.String),
        pl.col('percent_of_species_bacteria').map_elements(curr_format_float_number, return_dtype=pl.String),
        pl.col('number_of_species_archaea').map_elements(curr_format_int_number, return_dtype=pl.String),
        pl.col('percent_of_species_archaea').map_elements(curr_format_float_number, return_dtype=pl.String)
    )

    return fmt_copy_number_df
# end def
