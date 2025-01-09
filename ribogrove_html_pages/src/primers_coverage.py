
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

    # TODO: remove
    # df = df.rename(
    #     columns={
    #         'Taxon'                        : 'Phylum',
    #         'Number of genomes'            : 'num_genomes',
    #         '27F-1492R; Full gene (%)'     : '27F-1492R',
    #         '27F-338R; V1-V2 (%)'          : '27F-338R',
    #         '27F-534R; V1-V3 (%)'          : '27F-534R',
    #         '341F-785R; V3-V4 (%)'         : '341F-785R',
    #         '341F-944R; V3-V5 (%)'         : '341F-944R',
    #         '515F-806R; V4 (%)'            : '515F-806R',
    #         '515F-944R; V4-V5 (%)'         : '515F-944R',
    #         '515F-1100R; V4-V6 (%)'        : '515F-1100R',
    #         '784F-1100R; V5-V6 (%)'        : '784F-1100R',
    #         '784F-1193R; V5-V7 (%)'        : '784F-1193R',
    #         '939F-1193R; V6-V7 (%)'        : '939F-1193R',
    #         '939F-1378R; V6-V8 (%)'        : '939F-1378R',
    #         '1115F-1492R; V7-V9 (%)'       : '1115F-1492R',
    #         'SSU1ArF-SSU520R; V1-V4 (%)'   : 'SSU1ArF-SSU520R',
    #         '340f-806rB; V3-V4 (%)'        : '340f-806rB',
    #         '349f-519r; V3-V4 (%)'         : '349f-519r',
    #         '515fB-806rB; V4 (%)'          : '515fB-806rB',
    #         'Parch519f-Arch915r; V4-V5 (%)': 'Parch519f-Arch915r',
    #         '1106F-Ar1378R; V7-V8 (%)'     : '1106F-Ar1378R',
    #     }
    # )
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
    # TODO: remove
    # ordered_primer_names = [
    #     # Bacteria
    #     # Full
    #     '27F-1492R',
    #     # V1-V2
    #     '27F-338R',
    #     # V1-V3
    #     '27F-534R',
    #     # V3-V4
    #     '341F-785R',
    #     # V3-V5
    #     '341F-944R',
    #     # V4
    #     '515F-806R',
    #     # V4-V5
    #     '515F-944R',
    #     # V4-V6
    #     '515F-1100R',
    #     # V5-V6
    #     '784F-1100R',
    #     # V5-V7
    #     '784F-1193R',
    #     # V6-V7
    #     '939F-1193R',
    #     # V6-V8
    #     '939F-1378R',
    #     # V7-V9
    #     '1115F-1492R',
    #     # Archaea
    #     # V1-V4
    #     'SSU1ArF-SSU520R',
    #     # V3-V4
    #     '340f-806rB',
    #     # V3-V4
    #     '349f-519r',
    #     # V4
    #     '515fB-806rB',
    #     # V4-V5
    #     'Parch519f-Arch915r',
    #     # V7-V8
    #     '1106F-Ar1378R',
    # ]

    fmt_primer_df = primer_df[
        ['Domain', 'Phylum', 'num_genomes'] + ordered_primer_names
    ].sort_values(by='num_genomes', ascending=False)

    return fmt_primer_df.copy()
# end def


def _format_phylum_name(phylum_name):
    return phylum_name.replace('Candidatus ', 'Ca. ')
# end def

