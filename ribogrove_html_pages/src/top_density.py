
from functools import partial

import polars as pl

from src.formatting import format_int_number, format_float_number


def make_top_density_df(density_df,
                        source_genomes_df,
                        gene_stats_df,
                        top_type='highest',
                        top_num=10) -> pl.DataFrame:

    gcn_df = gene_stats_df.group_by('asm_acc').agg(
        pl.col('seqID').n_unique().alias('16S_rRNA_gcn')
    )

    density_df = density_df.join(
        gcn_df,
        on='asm_acc',
        how='left'
    ).join(
        source_genomes_df.select(pl.col('asm_acc', 'genome_size')),
        on='asm_acc',
        how='left'
    ).join(
        gene_stats_df.select(pl.col('asm_acc', 'strain_name', 'Domain')),
        on='asm_acc',
        how='left'
    ).unique()


    # Columns for the output dataframe
    out_columns = [
        'asm_acc',
        '16S_rRNA_density',
        '16S_rRNA_gcn',
        'genome_size',
        'strain_name',
        'Domain',
    ]

    sort_descending = top_type == 'highest'

    top_dfs = []

    # Do it for bacteria and for archaea
    for domain in ('Bacteria', 'Archaea'):
        domain_top_density_df = (
            density_df
            .filter(pl.col('Domain') == domain)
            .sort(by='16S_rRNA_density', descending=sort_descending)
            .head(top_num)
        )
        
        top_dfs.append(
            domain_top_density_df.select(pl.col(out_columns))
        )
    # end for

    top_df = pl.concat(top_dfs)

    print(top_df)

    return top_df
# end def

def _check_next_trait_the_same(genomes_at_trait_df, curr_top_i):
    curr_trait = genomes_at_trait_df.row(curr_top_i, named=True)['16S_rRNA_density']
    try:
        next_trait = genomes_at_trait_df.row(curr_top_i+1, named=True)['16S_rRNA_density']
    except pl.exceptions.OutOfBoundsError:
        # Catch case if number of rows is 1
        return False
    # end try
    return abs(curr_trait - next_trait) < 1e-9
# end def


def format_top_density_var_df(top_df, thousand_separator, decimal_separator):

    curr_format_int_number = partial(
        format_int_number,
        thousand_separator=thousand_separator
    )

    curr_format_float_number = partial(
        format_float_number,
        thousand_separator=thousand_separator,
        decimal_separator=decimal_separator,
        digits=3
    )

    fmt_top_df = top_df.with_columns(
        pl.col('16S_rRNA_density').map_elements(curr_format_float_number, return_dtype=pl.String),
        pl.col('16S_rRNA_gcn').map_elements(curr_format_int_number, return_dtype=pl.String),
        pl.col('genome_size').map_elements(curr_format_int_number, return_dtype=pl.String),
    )

    return fmt_top_df
# end def
