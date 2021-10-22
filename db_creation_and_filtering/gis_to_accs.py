#!/usr/bin/env python3

# The script takes RefSeq GI numbers from an input file (outputted by
#   the script `assembly2gi_numbers.py`) and finds corresponding RefSeq ACCESSION.VERSIONs and titles.
#   Output file is a TSV file of 3 columns (`gi_number`, `acc`, `title`).
# Requires Internet connection.
# For example, the script takes GI number 2075061612 (https://www.ncbi.nlm.nih.gov/nuccore/2075061612)
#   and finds corresponding ACCESSION.VERSION (NZ_CP063062.1)
#   and title ("Chlamydia suis strain 111 Ry chromosome, complete genome") for this record.

## Command line arguments
### Input files:
# 1. `-i / --gi-file` -- input TSV file mapping Assembly IDs to RefSeq GI numbers. Mandatory.

### Output files:
# 1. `-o / --outfile` -- output TSV file. Mandatory.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import time
import argparse

import pandas as pd

from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'


parser = argparse.ArgumentParser()

parser.add_argument(
    '-i',
    '--gi-file',
    help='TSV file (with header) with Assembly IDs and GI numbers separated by tabs',
    required=True
)

parser.add_argument(
    '-o',
    '--outfile',
    help='file mapping RefSeq GI numbers to corresponding ACCESSION.VERSION\'s and titles',
    required=True
)

args = parser.parse_args()

# For convenience
gi_fpath = os.path.realpath(args.gi_file)
outfpath = os.path.realpath(args.outfile)


# Check existance of input file
if not os.path.exists(gi_fpath):
    print(f'Error: file `{gi_fpath}` does not exist!')
    sys.exit(1)
# end if

if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

print(gi_fpath)
print()


# Read input
gi_df = pd.read_csv(
    gi_fpath,
    sep='\t'
)


chunk_size = 50
n_done_ids = 0

print('\r0/{}'.format(gi_df.shape[0]), end=' ')


# == Proceed ==
with open(outfpath, 'wt') as outfile:

    # Write header
    outfile.write('gi_number\tacc\ttitle\n')

    # Iterate over chunks of RefSeq IDs and get their "accession.version"s and titles
    for i in range(0, gi_df.shape[0], chunk_size):

        curr_gis = tuple(
            map(
                str,
                tuple(gi_df.iloc[i : i + chunk_size,]['gi_number'])
            )
        )

        # We will terminate if 3 errors occur in a row
        # Request summary for RefSeq record
        error = True
        n_errors = 0
        while error:
            try:
                handle = Entrez.esummary(
                    db='nuccore',
                    id=','.join(curr_gis)
                )
                records = Entrez.read(handle)
                handle.close()
            except OSError:
                n_errors += 1
                if n_errors == 3:
                    raise OSError("Cannot request esummary")
                # end if
            else:
                error = False
            # end try
        # end while

        # Rwite output
        for rec in records:
            outfile.write(f'{rec["Id"]}\t{rec["AccessionVersion"]}\t{rec["Title"]}\n')
        # end for

        # Wait a bit: we don't want NCBI to ban us :)
        time.sleep(0.4)

        n_done_ids += chunk_size
        print('\r{}/{}'.format(n_done_ids, gi_df.shape[0]), end=' ')
    # end for
# end with

print('\r{}/{}   '.format(min(n_done_ids, gi_df.shape[0]), gi_df.shape[0]))

print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
