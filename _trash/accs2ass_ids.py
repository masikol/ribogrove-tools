#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# The script takes IDs of the NCBI Assembly database from an input file and finds
#   corresponding RefSeq GI numbers using elink utility (https://www.ncbi.nlm.nih.gov/books/NBK25497/).
#   The script writes result to the output TSV file of 2 columns (`ass_id`, `gi_numbers`),
#   thus mapping Assembly IDs to RefSeq GI numbers.
# Requires Internet connection.
# For example, the script takes Assembly ID 10601591 (https://www.ncbi.nlm.nih.gov/assembly/10601591/)
#   and finds corresponding RefSeq record [2075061612](https://www.ncbi.nlm.nih.gov/nuccore/2075061612).

## Command line arguments

### Input files:
# 1. `-i / --accs-file` -- input file of Assembly IDs (one per line). Mandatory.

###  Output files:
# 1. `-o / --outfile` -- output TSV file mapping Assembly IDs to RefSeq GI numbers. Mandatory.

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import time
import argparse
import http.client


from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'

# == Parse arguments ==

parser = argparse.ArgumentParser()

parser.add_argument(
    '-i',
    '--accs-file',
    help='file with Assembly IDs, one per line',
    required=True
)

parser.add_argument(
    '-o',
    '--outfile',
    help='file mapping Assembly IDs to RefSeq GI numbers',
    required=True
)

args = parser.parse_args()


# For convenience
accs_fpath = os.path.realpath(args.accs_file)
outfpath = os.path.realpath(args.outfile)


# Check existance of input file
if not os.path.exists(accs_fpath):
    print(f'Error: file `{accs_fpath}` does not exist!')
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


# Read assembly IDs
accs = tuple(
    map(
        str.strip,
        open(accs_fpath).readlines()
    )
)

# == Proceed ==

with open(outfpath, 'wt') as outfile:

    # Write header
    # outfile.write('acc\tass_id\n')

    # Iterate over assembly IDs
    for i, acc in enumerate(accs):

        print(f'\rDoing {i+1}/{len(accs)}: {acc}', end=' ')

        # We will terminate if 3 errors occur in a row
        n_errors = 0

        # Request RefSeq GI number
        while n_errors < 3:
            try:
                handle = Entrez.elink(
                    dbfrom='nuccore',
                    db='assembly',
                    id=acc
                )
                record = Entrez.read(handle)
                handle.close()
                break
            except (IOError, RuntimeError, http.client.HTTPException) as err:
                print('\n' + str(err))
                n_errors += 1
            # end try

            if len(record[0]['ERROR']) != 0:
                print('\n' + '\n'.join(record[0]['ERROR']))
                n_errors += 1
            # end if
        # end while

        if n_errors >= 3:
            print('\n3 errors!')
            continue
        # end if

        # Extract RefSeq GI number
        try:
            ass_ids = [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']]
        except IndexError:
            # print(f'\nError on {acc}: {record}')
            for ass_id in ass_ids:
                outfile.write(f'{acc}\tNA\n')
            # end for
        else:
            # Write output
            for ass_id in ass_ids:
                outfile.write(f'{acc}\t{ass_id}\n')
            # end for
        # end try

        # Wait a bit: we don't want NCBI to ban us :)
        time.sleep(0.4)
    # end for
# end with

print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
