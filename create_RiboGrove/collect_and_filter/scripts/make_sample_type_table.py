#!/usr/bin/env python3

import os
from src.rg_tools_time import get_time

print(
    '\n|=== {} STARTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)


# == Parse arguments ==
import argparse

parser = argparse.ArgumentParser()

# Input
parser.add_argument(
    '-i',
    '--in-asm-data',
    help='The file `tmp_asm_data.json` produced by `datasets summary genome accession` program call',
    required=True
)

parser.add_argument(
    '-o',
    '--outfile',
    help='Output TSV table path',
    required=True
)


args = parser.parse_args()

infpath = os.path.abspath(args.in_asm_data)
outfpath = os.path.abspath(args.outfile)


if not os.path.isfile(infpath):
    print(f'Error: file `{asm_sum_fpath}` does not exist')
    sys.exit(1)
# end if


# == Import them now ==

import sys
import json

# == Proceed ==

with open(infpath, 'rt') as in_handle:
    asm_dict = json.load(in_handle)
# end with
sep = '\t'

with open(outfpath, 'wt') as out_handle:

    out_handle.write('{}\n'.format(sep.join([
        'asm_acc',
        'sample_type'
    ])))

    for report in asm_dict['reports']:
        asm_acc = report['accession']

        non_available = False

        try:
            attributes = report['assembly_info']['biosample']['attributes']
        except KeyError:
            sample_type='NA'
            non_available = True
        # end try

        if not non_available:
            sample_type = 'NA'
            try:
                for attr_dict in attributes:
                    if attr_dict['name'] == 'sample_type':
                        sample_type = attr_dict['value']
                    # end if
                # end for
            # end try
            except KeyError:
                pass
            # end try
        # end if

        out_handle.write('{}\n'.format(sep.join([
            asm_acc,
            sample_type
        ])))
    # end for
# end with


print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)

sys.exit(0)
