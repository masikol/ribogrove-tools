from functools import partial

import pandas as pd

from src.formatting import format_int_number, format_float_number

def make_ribogrove_top_intragenomic_var_df(entropy_summary_df, gene_stats_df, top_num=10):

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
    series_nunique = lambda x: x.nunique()
    by_genome_copy_number_df = gene_stats_df.groupby('asm_acc', as_index=False) \
        .agg({'seqID': series_nunique}) \
        .rename(columns={'seqID': 'copy_number'}) \
        .merge(
            gene_stats_df[['asm_acc', 'strain_name', 'Domain']].drop_duplicates(),
            on='asm_acc',
            how='left'
        )

    # Map Assembly IDs to domain names
    entropy_summary_df = entropy_summary_df.merge(
        by_genome_copy_number_df,
        on='asm_acc',
        how='left'
    )

    # Create an output dataframe
    top_df = pd.DataFrame({colname: [] for colname in out_columns})

    # Do it for bacteria and for archaea
    for domain in ('Bacteria', 'Archaea'):

        # Create a dataframe of maximum sum of entropy for each genome
        domain_entropy_df = entropy_summary_df[
            entropy_summary_df['Domain'] == domain
        ].sort_values(by='sum_entropy', ascending=False) \
            .reset_index()

        if domain_entropy_df.shape[0] == 0:
            return pd.DataFrame(
                {
                    'asm_acc': list(),
                    'sum_entropy': list(),
                    'mean_entropy': list(),
                    'num_var_cols': list(),
                    'copy_number': list(),
                    'strain_name': list(),
                    'Domain': list(),
                }
            )
        # end if

        # We will stop if we reach `top_num` genomes
        #   and if the next (`top_num`+1)th one has the same sum of entropy as `top_num`th one
        #   we wil add (`top_num`+1)th genome too
        top_genome_counter = 0
        next_sum_entropy_the_same = abs(
            domain_entropy_df.loc[top_genome_counter, 'sum_entropy'] \
            - domain_entropy_df.loc[top_genome_counter+1, 'sum_entropy']
        ) < 1e-9

        while top_genome_counter < top_num or next_sum_entropy_the_same:

            series_to_append = pd.Series(
                {
                    'asm_acc': domain_entropy_df.loc[top_genome_counter, 'asm_acc'],
                    'sum_entropy': domain_entropy_df.loc[top_genome_counter, 'sum_entropy'],
                    'mean_entropy': domain_entropy_df.loc[top_genome_counter, 'mean_entropy'],
                    'num_var_cols': domain_entropy_df.loc[top_genome_counter, 'num_var_cols'],
                    'copy_number': domain_entropy_df.loc[top_genome_counter, 'copy_number'],
                    'strain_name': domain_entropy_df.loc[top_genome_counter, 'strain_name'],
                    'Domain': domain_entropy_df.loc[top_genome_counter, 'Domain'],
                }
            )

            # Append selected rows to the output dataframe
            top_df = top_df.append(series_to_append, ignore_index=True)

            if top_genome_counter == domain_entropy_df.shape[0] - 1:
                break
            # end if
            # Update contition variables
            next_sum_entropy_the_same = abs(
                domain_entropy_df.loc[top_genome_counter, 'sum_entropy'] \
                - domain_entropy_df.loc[top_genome_counter+1, 'sum_entropy']
            ) < 1e-9
            top_genome_counter += 1
        # end while
    # end for

    top_df = top_df.reset_index()
    top_df['num_var_cols'] = top_df['num_var_cols'].map(int)

    print(top_df)

    return top_df
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

    fmt_top_df = top_df.copy()
    fmt_top_df['sum_entropy'] = fmt_top_df['sum_entropy'] \
        .map(curr_format_float_number)
    fmt_top_df['mean_entropy'] = fmt_top_df['mean_entropy'] \
        .map(curr_format_float_number)
    fmt_top_df['num_var_cols'] = fmt_top_df['num_var_cols'] \
        .map(curr_format_int_number)
    fmt_top_df['copy_number'] = fmt_top_df['copy_number'] \
        .map(curr_format_int_number)

    return fmt_top_df
# end def
