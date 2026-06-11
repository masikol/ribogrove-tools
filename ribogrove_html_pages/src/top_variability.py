
from functools import partial

import polars as pl

from src.formatting import format_int_number, format_float_number


_IGNORE_ASM_ACCS = {
    'GCF_019974355.1', # just a large insertion in GCF_019974355.1:NZ_AP024929.1:249100-251537:minus. other positions are identical
}


def make_ribogrove_top_intragenomic_var_df(entropy_summary_df, gene_stats_df, top_num=10):

    # TODO: remove
    # entropy_tmp_df = entropy_summary_df.query('not asm_acc in @_IGNORE_ASM_ACCS')
    entropy_tmp_df = entropy_summary_df.filter(
        ~pl.col('asm_acc').is_in(_IGNORE_ASM_ACCS)
    )

    # Columns for the output dataframe
    out_columns = [
        'asm_acc',
        'sum_entropy',
        'mean_entropy',
        'num_var_cols',
        'copy_number',
        'strain_name',
        'Domain',
    ]

    # Count copy numbers
    # TODO: remove pd
    # series_nunique = lambda x: x.nunique()
    # gcn_df = gene_stats_df.groupby('asm_acc', as_index=False) \
    #     .agg({'seqID': series_nunique}) \
    #     .rename(columns={'seqID': 'copy_number'}) \
    #     .merge(
    #         gene_stats_df[['asm_acc', 'strain_name', 'Domain']].drop_duplicates(),
    #         on='asm_acc',
    #         how='left'
    #     )
    gcn_df = gene_stats_df.group_by('asm_acc').agg(
        pl.col('seqID').n_unique().alias('copy_number')
    ).join(
        gene_stats_df.select(pl.col('asm_acc', 'strain_name', 'Domain')).unique(),
        on='asm_acc',
        how='left'
    )

    # Map Assembly IDs to domain names
    # TODO: remove pd
    # entropy_tmp_df = entropy_tmp_df.merge(
    #     gcn_df,
    #     on='asm_acc',
    #     how='left'
    # )
    entropy_tmp_df = entropy_tmp_df.join(
        gcn_df,
        on='asm_acc',
        how='left'
    )

    # Create an output dataframe
    # TODO: remove pd
    # top_df = pd.DataFrame({colname: [] for colname in out_columns})
    top_rows = []

    # Do it for bacteria and for archaea
    for domain in ('Bacteria', 'Archaea'):

        # Create a dataframe of maximum sum of entropy for each genome
        # TODO: remove pd
        # domain_entropy_df = entropy_tmp_df[
        #     entropy_tmp_df['Domain'] == domain
        # ].sort_values(by='sum_entropy', ascending=False) \
        #     .reset_index()
        domain_entropy_df = entropy_tmp_df \
            .filter(pl.col('Domain') == domain) \
            .sort(by='sum_entropy', descending=True)

        if domain_entropy_df.height == 0:
            return pl.DataFrame(
                {colname: list() for colname in out_columns}
            )
        # end if

        # We will stop if we reach `top_num` genomes
        #   and if the next (`top_num`+1)th one has the same sum of entropy as `top_num`th one
        #   we wil add (`top_num`+1)th genome too
        top_i = 0
        next_sum_entropy_the_same = _check_next_trait_the_same(domain_entropy_df, top_i)

        while top_i < top_num or next_sum_entropy_the_same:

            row = domain_entropy_df.row(top_i, named=True)

            top_rows.append({
                'asm_acc':      row['asm_acc'],
                'sum_entropy':  row['sum_entropy'],
                'mean_entropy': row['mean_entropy'],
                'num_var_cols': row['num_var_cols'],
                'copy_number':  row['copy_number'],
                'strain_name':  row['strain_name'],
                'Domain':       row['Domain'],
            })

            if top_i == domain_entropy_df.height - 1:
                break
            # end if

            # Update contition variables
            next_sum_entropy_the_same = _check_next_trait_the_same(domain_entropy_df, top_i)
            top_i += 1
        # end while
    # end for

    # Build final DataFrame
    if len(top_rows) != 0:
        top_df = pl.DataFrame(top_rows, orient='row')
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

def _check_next_trait_the_same(genomes_at_trait_df, curr_top_i):
    curr_trait = genomes_at_trait_df.row(curr_top_i, named=True)['sum_entropy']
    try:
        next_trait = genomes_at_trait_df.row(curr_top_i+1, named=True)['sum_entropy']
    except pl.exceptions.OutOfBoundsError:
        # Catch case if number of rows is 1
        return False
    # end try
    return abs(curr_trait - next_trait) < 1e-9
# end def


def format_top_intragenomic_var_df(top_df, thousand_separator, decimal_separator):

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

    fmt_top_df = top_df.with_columns(
        pl.col('sum_entropy').map_elements(curr_format_float_number, return_dtype=pl.String),
        pl.col('mean_entropy').map_elements(curr_format_float_number, return_dtype=pl.String),
        pl.col('num_var_cols').map_elements(curr_format_int_number, return_dtype=pl.String),
        pl.col('copy_number').map_elements(curr_format_int_number, return_dtype=pl.String)
    )

    return fmt_top_df
# end def
