
import os
import json
from collections import OrderedDict
from functools import partial, reduce

import pandas as pd

from src.formatting import format_int_number, format_float_number


def make_ribogrove_primer_coverage_df(input_fpath):
    df = pd.read_csv(input_fpath, sep='\t')
    df = df[df['Rank'] == 'Phylum']

    df = df.drop(['Rank'], axis=1)

    bacterial_primer_pairs, archaeal_primer_pairs = parse_primer_pairs()

    renamed_columns = {
        'Taxon'                        : 'Phylum',
        'Number of genomes'            : 'num_genomes',
    }
    for domain_primer_pairs in (bacterial_primer_pairs, archaeal_primer_pairs):
        d_rename_columns = {
            '{}; {} (%)'.format(name_pair, v_region_name): name_pair
                for name_pair, v_region_name in domain_primer_pairs.items()
        }
        renamed_columns.update(d_rename_columns)
    # end for

    df = df.rename(
        columns=renamed_columns
    )

    print(df)

    return df[
        [
            'Domain',
            'Phylum',
            'num_genomes',
        ] \
        + list(bacterial_primer_pairs.keys()) \
        + list(archaeal_primer_pairs.keys())
    ]
# end def


def parse_primer_pairs():
    # TODO: deduplicate code
    primers_pairs_fpath = os.path.join(
        os.path.dirname(os.path.dirname(os.path.dirname(__file__))),
        'create_RiboGrove', 'collect_and_filter', 'scripts', 'data', 'primers', 'primer_pairs.json'
    )
    with open(primers_pairs_fpath, 'rt') as infile:
        primer_pairs = json.load(infile)
    # end with

    bacterial_primer_pairs = transform_primer_pair_dict(primer_pairs['Bacteria'])
    archaeal_primer_pairs  = transform_primer_pair_dict(primer_pairs['Archaea'])

    return bacterial_primer_pairs, archaeal_primer_pairs
# end def

def transform_primer_pair_dict(primer_pairs):
    all_primer_pair_dict = OrderedDict()
    for nameF, nameR, v_region_name in primer_pairs:
        all_primer_pair_dict['{}-{}'.format(nameF, nameR)] = v_region_name
    # end for
    return all_primer_pair_dict
# end def


def format_primer_coverage_df(primer_coverage_df,
                              thousand_separator,
                              decimal_separator):

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

    fmt_primer_df = _sort_rows_and_columns(primer_coverage_df)

    fmt_primer_df['Phylum'] = fmt_primer_df['Phylum'] \
        .map(_format_phylum_name)

    fmt_primer_df['num_genomes'] = fmt_primer_df['num_genomes'] \
        .map(curr_format_int_number)

    bacterial_primer_pairs, archaeal_primer_pairs = parse_primer_pairs()
    primer_keys = list(bacterial_primer_pairs.keys()) \
                  + list(archaeal_primer_pairs.keys())

    for primer_key in primer_keys:
        fmt_primer_df[primer_key] = fmt_primer_df[primer_key] \
            .map(curr_format_float_number)
    # end for


    return fmt_primer_df
# end def


def _sort_rows_and_columns(primer_df):
    bacterial_primer_pairs, archaeal_primer_pairs = parse_primer_pairs()
    ordered_primer_names = list(bacterial_primer_pairs.keys()) \
                         + list(archaeal_primer_pairs.keys())
    fmt_primer_df = primer_df[
        ['Domain', 'Phylum', 'num_genomes'] + ordered_primer_names
    ].sort_values(by='num_genomes', ascending=False)

    return fmt_primer_df.copy()
# end def


def _format_phylum_name(phylum_name):
    return phylum_name.replace('Candidatus ', 'Ca. ')
# end def

