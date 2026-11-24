
from functools import partial

import polars as pl

from src.util import get_top_rows
from src.formatting import format_int_number


def make_ribogrove_top_longest_df(gene_stats_df: pl.DataFrame, top_num: int = 10):
    out_columns = ['asm_acc', 'len', 'seqID', 'strain_name', 'Domain']
    top_rows = []

    for domain in ('Bacteria', 'Archaea'):
        # Filter by domain and get max length per asm_acc
        domain_df = gene_stats_df.filter(pl.col('Domain') == domain)
        
        max_len_domain_df = (
            domain_df
            .group_by('asm_acc')
            .agg(pl.col('len').max().alias('len'))
            .sort('len', descending=True)
        )
            
        # Get unique lengths for tie handling
        unique_lengths = max_len_domain_df['len'].unique().sort(descending=True)

        top_rows.extend(
            get_top_rows(
                domain_df,
                max_len_domain_df,
                unique_lengths,
                trait_col_name='len',
                top_num=top_num
            )
        )
    # end for

    # Build final DataFrame
    if len(top_rows) != 0:
        top_df = pl.DataFrame(top_rows, orient='row')
        top_df = top_df.with_columns(pl.col('len').cast(pl.UInt64))
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


def format_longest_genes_df(top_df, thousand_separator, decimal_separator):

    curr_format_int_number = partial(
        format_int_number,
        thousand_separator=thousand_separator
    )

    fmt_top_df = top_df.with_columns(
        pl.col('len').map_elements(curr_format_int_number, return_dtype=pl.String)
    )

    return fmt_top_df
# end def
