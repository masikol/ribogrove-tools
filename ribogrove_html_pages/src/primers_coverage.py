
import os
import json
from functools import partial
from collections import OrderedDict

import polars as pl

from src.formatting import format_int_number, format_float_number


def make_ribogrove_primer_coverage_df(input_fpath):
    df = (
        pl.scan_csv(input_fpath, separator='\t')
        .filter(pl.col('Rank') == 'Phylum')
        .drop('Rank')
        .collect()
    )

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

    df = df.rename(renamed_columns)

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
    return df.select(pl.col(
        ['Domain', 'Phylum', 'num_genomes',] \
        + list(bacterial_primer_pairs.keys()) \
        + list(archaeal_primer_pairs.keys())
    ))
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


    fmt_primer_df = fmt_primer_df.with_columns(
        pl.col("Phylum").str.replace_all('Candidatus ','Ca. ',literal=True),
        pl.col('num_genomes').map_elements(curr_format_int_number, return_dtype=pl.String),
    )

    bacterial_primer_pairs, archaeal_primer_pairs = parse_primer_pairs()
    primer_keys = list(bacterial_primer_pairs.keys()) \
                  + list(archaeal_primer_pairs.keys())

    for primer_key in primer_keys:
        fmt_primer_df = fmt_primer_df.with_columns(
            pl.col(primer_key).map_elements(curr_format_float_number, return_dtype=pl.String)
        )
    # end for

    return fmt_primer_df
# end def


def _sort_rows_and_columns(primer_df):
    bacterial_primer_pairs, archaeal_primer_pairs = parse_primer_pairs()
    ordered_primer_names = list(bacterial_primer_pairs.keys()) \
                         + list(archaeal_primer_pairs.keys())
    fmt_primer_df = primer_df.select(pl.col(
        ['Domain', 'Phylum', 'num_genomes'] + ordered_primer_names
    )).sort(by='num_genomes', descending=True)

    return fmt_primer_df
# end def


def _format_phylum_name(phylum_name):
    return phylum_name.replace('Candidatus ', 'Ca. ')
# end def
