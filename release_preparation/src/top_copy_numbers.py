
from functools import partial

import pandas as pd

from src.formatting import format_int_number


def make_ribogrove_top_copy_numbers_df(gene_stats_df, top_num=10):

    series_nunique = lambda x: x.nunique()

    # Columns for the output dataframe
    out_columns = ['ass_id', 'copy_number', 'title', 'domain']

    by_genome_copy_number_df = gene_stats_df.groupby('ass_id', as_index=False) \
        .agg({'seqID': series_nunique}) \
        .rename(columns={'seqID': 'copy_number'}) \
        .merge(
            gene_stats_df[['ass_id', 'title', 'domain']].drop_duplicates(),
            on='ass_id',
            how='left'
        )

    # Create an output dataframe
    top_df = pd.DataFrame({colname: [] for colname in out_columns})

    # Do it for bacteria and for archaea
    for domain in ('Bacteria', 'Archaea'):

        # Create a dataframe of maximum gene copy number for each genome
        domain_copy_number_df = by_genome_copy_number_df[
            by_genome_copy_number_df['domain'] == domain
        ].sort_values(by='copy_number', ascending=False) \
            .reset_index()

        # We will stop if we reach `top_num` genomes
        #   and if the next (`top_num`+1)th one has the same maximum gene copy number as `top_num`th one
        #   we wil add (`top_num`+1)th genome too
        top_genome_counter = 0
        next_copy_number_the_same = domain_copy_number_df.loc[top_genome_counter, 'copy_number'] \
                            == domain_copy_number_df.loc[top_genome_counter+1, 'copy_number']

        while top_genome_counter < top_num or next_copy_number_the_same:

            series_to_append = pd.Series(
                {
                    'ass_id': domain_copy_number_df.loc[top_genome_counter, 'ass_id'],
                    'copy_number': domain_copy_number_df.loc[top_genome_counter, 'copy_number'],
                    'title': domain_copy_number_df.loc[top_genome_counter, 'title'],
                    'domain': domain_copy_number_df.loc[top_genome_counter, 'domain'],
                }
            )

            # Append selected rows to the output dataframe
            top_df = top_df.append(series_to_append, ignore_index=True)

            # Update contition variables
            next_copy_number_the_same = domain_copy_number_df.loc[top_genome_counter, 'copy_number'] \
                                == domain_copy_number_df.loc[top_genome_counter+1, 'copy_number']
            top_genome_counter += 1
        # end while
    # end for

    top_df = top_df.reset_index()

    top_df['ass_id'] = top_df['ass_id'].map(int)
    top_df['copy_number'] = top_df['copy_number'].map(int)

    print(top_df)

    return top_df
# end def make_ribogrove_top_copy_numbers_df


def format_top_copy_numbers_df(top_df, thousand_separator, decimal_separator):

    curr_format_int_number = partial(
        format_int_number,
        thousand_separator=thousand_separator
    )

    fmt_top_df = top_df.copy()
    fmt_top_df['copy_number'] = fmt_top_df['copy_number'] \
        .map(curr_format_int_number)

    return fmt_top_df
# end def format_top_copy_numbers_df
