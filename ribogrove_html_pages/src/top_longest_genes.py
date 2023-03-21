
from functools import partial

import pandas as pd

from src.formatting import format_int_number


def make_ribogrove_top_longest_df(gene_stats_df, top_num=10):
    # Columns for the output dataframe
    out_columns = ['asm_acc', 'len', 'seqID', 'strain_name', 'Domain']

    # Create an output dataframe
    top_df = pd.DataFrame({colname: [] for colname in out_columns})

    # Do it for bacteria and for archaea
    for domain in ('Bacteria', 'Archaea'):

        # Create a dataframe of maximum gene length for each genome
        max_len_domain_df = gene_stats_df[
            gene_stats_df['Domain'] == domain
        ].groupby('asm_acc', as_index=False).agg({'len': 'max'}) \
            .sort_values(by='len', ascending=False) \
            .reset_index()

        # We will stop if we reach `top_num` genomes
        #   and if the next (`top_num`+1)th one has the same maximum gene length as `top_num`th one
        #   we wil add (`top_num`+1)th genome too
        top_genome_counter = 0
        next_len_the_same = max_len_domain_df.loc[top_genome_counter, 'len'] \
                            == max_len_domain_df.loc[top_genome_counter+1, 'len']

        while top_genome_counter < top_num or next_len_the_same:

            # Get current Assembly ID and it's maxumum gene length
            curr_asm_acc = max_len_domain_df.loc[top_genome_counter, 'asm_acc']
            curr_max_len = max_len_domain_df.loc[top_genome_counter, 'len']

            # Get all row from `curr_genome_df` of the current Assembly
            #   and of maximum gene length
            curr_genome_df = gene_stats_df[
                (gene_stats_df['asm_acc'] == curr_asm_acc) \
                & (gene_stats_df['len'] == curr_max_len)
            ][out_columns]

            series_to_append = pd.Series(
                {
                    'asm_acc': list(curr_genome_df['asm_acc'])[0],
                    'len': list(curr_genome_df['len'])[0],
                    'seqID': list(curr_genome_df['seqID']),
                    'strain_name': list(curr_genome_df['strain_name'])[0],
                    'Domain': list(curr_genome_df['Domain'])[0],
                }
            )

            # Append selected rows to the output dataframe
            top_df = pd.concat(
                [
                    top_df,
                    series_to_append.to_frame().T,
                ],
                ignore_index=True
            )

            if top_genome_counter == max_len_domain_df.shape[0] - 1:
                break
            # end if
            next_len_the_same = max_len_domain_df.loc[top_genome_counter, 'len'] \
                                == max_len_domain_df.loc[top_genome_counter+1, 'len']
            # Update contition variables
            top_genome_counter += 1
        # end while
    # end for

    top_df = top_df.reset_index()
    top_df['len'] = top_df['len'].map(int)

    print(top_df)

    return top_df
# end def


def format_longest_genes_df(top_df, thousand_separator, decimal_separator):

    curr_format_int_number = partial(
        format_int_number,
        thousand_separator=thousand_separator
    )

    fmt_top_df = top_df.copy()
    fmt_top_df['len'] = fmt_top_df['len'] \
        .map(curr_format_int_number)

    return fmt_top_df
# end def
