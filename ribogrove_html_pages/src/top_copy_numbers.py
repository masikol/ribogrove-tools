
from functools import partial

import polars as pl

from src.util import get_top_rows
from src.formatting import format_int_number


def make_ribogrove_top_copy_numbers_df(gene_stats_df: pl.DataFrame, top_num: int = 10):
    out_columns = ['asm_acc', 'copy_number', 'seqID', 'strain_name', 'Domain']
    top_rows = []


    for domain in ('Bacteria', 'Archaea'):
        # Filter by domain and get max copy number per asm_acc
        domain_df = gene_stats_df.filter(pl.col('Domain') == domain)

        gcn_df = domain_df.group_by('asm_acc').agg(
            pl.col('seqID').n_unique().alias('copy_number')
        )

        gcn_domain_df = (
            domain_df
            .join(gcn_df, on='asm_acc', how='left')
        )


        # Get unique numbers for tie handling
        unique_copy_numbers = (
            gcn_domain_df['copy_number']
            .unique()
            .sort(descending=True)
            .to_list()
        )

        top_rows.extend(
            get_top_rows(
                gcn_domain_df,
                gcn_df,
                unique_copy_numbers,
                trait_col_name='copy_number',
                top_num=top_num
            )
        )
    # end for

    # Build final DataFrame
    if len(top_rows) != 0:
        top_df = pl.DataFrame(top_rows, orient='row')
        top_df = top_df.with_columns(pl.col('copy_number').cast(pl.UInt64))
    else:
        top_df = pl.DataFrame(
            schema={
                col: pl.String() for col in out_columns
            }
        )
    # end if
    
    print(top_df)
    return top_df
# end def


def format_top_copy_numbers_df(top_df, thousand_separator, decimal_separator):

    curr_format_int_number = partial(
        format_int_number,
        thousand_separator=thousand_separator
    )

    fmt_top_df = top_df.with_columns(
        pl.col('copy_number').map_elements(curr_format_int_number, return_dtype=pl.String)
    )

    return fmt_top_df
# end def
