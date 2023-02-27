
from functools import partial

import numpy as np
import pandas as pd

from src.formatting import format_float_number

def _init_len_dict():

    ribogrove_len_dict = {
        'min': {
            'Bacteria': None,
            'Archaea':  None,
        },
        '25perc': {
            'Bacteria': None,
            'Archaea':  None,
        },
        'median': {
            'Bacteria': None,
            'Archaea':  None,
        },
        '75perc': {
            'Bacteria': None,
            'Archaea':  None,
        },
        'mean': {
            'Bacteria': None,
            'Archaea':  None,
        },
        'modes': {
            'Bacteria': None,
            'Archaea':  None,
        },
        'max': {
            'Bacteria': None,
            'Archaea':  None,
        },
        'std': {
            'Bacteria': None,
            'Archaea':  None,
        },
    }

    return ribogrove_len_dict
# end def _init_len_dict


def make_ribogrove_len_dict(gene_stats_df):
    ribogrove_len_dict = _init_len_dict()

    bacteria_df = gene_stats_df[
        gene_stats_df['domain'] == 'Bacteria'
    ]
    archaea_df = gene_stats_df[
        gene_stats_df['domain'] == 'Archaea'
    ]

    # Calculate minimum lengths
    ribogrove_len_dict['min']['Bacteria'] = bacteria_df['len'].min()
    ribogrove_len_dict['min']['Archaea'] = archaea_df['len'].min()

    # Calculate maximum lengths
    ribogrove_len_dict['max']['Bacteria'] = bacteria_df['len'].max()
    ribogrove_len_dict['max']['Archaea'] = archaea_df['len'].max()

    # Make normalized (by species) dataframe
    bacteria_normalized_df = bacteria_df.groupby('species', as_index=False) \
        .agg({'len': 'median'})
    archaea_normalized_df = archaea_df.groupby('species', as_index=False) \
        .agg({'len': 'median'})

    # 25-th percentile
    ribogrove_len_dict['25perc']['Bacteria'] = np.percentile(bacteria_normalized_df['len'], 25)
    ribogrove_len_dict['25perc']['Archaea'] = np.percentile(archaea_normalized_df['len'], 25)

    # Median
    ribogrove_len_dict['median']['Bacteria'] = bacteria_normalized_df['len'].median()
    ribogrove_len_dict['median']['Archaea'] = archaea_normalized_df['len'].median()

    # 75-th percentile
    ribogrove_len_dict['75perc']['Bacteria'] = np.percentile(bacteria_normalized_df['len'], 75)
    ribogrove_len_dict['75perc']['Archaea'] = np.percentile(archaea_normalized_df['len'], 75)

    # Mean
    ribogrove_len_dict['mean']['Bacteria'] = bacteria_normalized_df['len'].mean()
    ribogrove_len_dict['mean']['Archaea'] = archaea_normalized_df['len'].mean()

    # Mode. There can bu multiple modes for a distribution
    bacteria_modes = list(bacteria_normalized_df['len'].mode())
    ribogrove_len_dict['modes']['Bacteria'] = bacteria_modes
    archaea_modes = list(archaea_normalized_df['len'].mode())
    ribogrove_len_dict['modes']['Archaea'] = archaea_modes

    # Standard deviation
    ribogrove_len_dict['std']['Bacteria'] = bacteria_normalized_df['len'].std()
    ribogrove_len_dict['std']['Archaea'] = archaea_normalized_df['len'].std()

    for k, v in ribogrove_len_dict.items():
        print(f'{k}: {v}')
    # end for

    return ribogrove_len_dict
# end def make_ribogrove_len_dict


def format_len_dict(ribogrove_len_dict, thousand_separator, decimal_separator):
    fmt_ribogrove_len_dict = dict()

    curr_format_float = partial(
        format_float_number,
        thousand_separator=thousand_separator,
        decimal_separator=decimal_separator,
        digits=2
    )

    for metric_name in ribogrove_len_dict.keys():
        fmt_ribogrove_len_dict[metric_name] = dict()
        for column_name, metric_value in ribogrove_len_dict[metric_name].items():
            if metric_name != 'modes':
                fmt_ribogrove_len_dict[metric_name][column_name] = curr_format_float(metric_value)
            else:
                # Format modes. There can be multiple modes
                fmt_ribogrove_len_dict[metric_name][column_name] = '<br />'.join(
                    map(
                        curr_format_float,
                        metric_value
                    )
                )
            # end if
        # end for
    # end for
    return fmt_ribogrove_len_dict
# end def format_len_dict
