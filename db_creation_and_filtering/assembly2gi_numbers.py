#!/usr/bin/env python3

# The script takes IDs of the NCBI Assembly database from an input file and finds
#   corresponding RefSeq GI numbers using elink utility (https://www.ncbi.nlm.nih.gov/books/NBK25497/).
#   The script writes result to the output TSV file of 2 columns (`ass_id`, `gi_numbers`),
#   thus mapping Assembly IDs to RefSeq GI numbers.
# Requires Internet connection.
# For example, the script takes Assembly ID 10601591 (https://www.ncbi.nlm.nih.gov/assembly/10601591/)
#   and finds corresponding RefSeq record [2075061612](https://www.ncbi.nlm.nih.gov/nuccore/2075061612).

## Command line arguments

### Input files:
# 1. `-i / --assm-id-file` -- input file of Assembly IDs (one per line). Mandatory.

###  Output files:
# 1. `-o / --outfile` -- output TSV file mapping Assembly IDs to RefSeq GI numbers. Mandatory.

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import time
import argparse
import http.client

import pandas as pd

from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'

# == Parse arguments ==

parser = argparse.ArgumentParser()

parser.add_argument(
    '-i',
    '--assm-id-file',
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
assm_id_fpath = os.path.realpath(args.assm_id_file)
outfpath = os.path.realpath(args.outfile)


# Check existance of input file
if not os.path.exists(assm_id_fpath):
    print(f'Error: file `{assm_id_fpath}` does not exist!')
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
ass_ids = tuple(
    map(
        str.strip,
        open(assm_id_fpath).readlines()
    )
)


# == Proceed ==

with open(outfpath, 'wt') as outfile:

    # Write header
    outfile.write('ass_id\tgi_number\n')

    # Iterate over assembly IDs
    for i, ass_id in enumerate(ass_ids):

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' ')

        # We will terminate if 3 errors occur in a row
        n_errors = 0

        # Request RefSeq GI number
        while n_errors < 3:
            try:
                handle = Entrez.elink(
                    dbfrom='assembly',
                    db='nuccore',
                    id=ass_id,
                    linkname='assembly_nuccore_refseq'
                )
                record = Entrez.read(handle)
                handle.close()
                break
            except (IOError, RuntimeError, http.client.IncompleteRead) as err:
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
            print(f'Please, handle Assembly ID {ass_id} manually.')
            print('Sleeping for 60 seconds')
            time.sleep(60)
            continue
        # end if

        # Extract RefSeq GI number
        try:
            gi_numbers = [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']]
        except IndexError:
            pass
        else:
            # Write output
            for gi_number in gi_numbers:
                outfile.write(f'{ass_id}\t{gi_number}\n')
            # end for
        # end try

        # Wait a bit: we don't want NCBI to ban us :)
        time.sleep(0.4)
    # end for
# end with

print('\nRemoving duplicated GI numbers...')

raw_df = pd.read_csv(outfpath, sep='\t')
drop_dupl_df = raw_df.drop_duplicates(subset=['gi_number'])

if drop_dupl_df.shape[0] != raw_df.shape[0]:
    print(
        'Found {} duplicated GI numbers' \
            .format(raw_df.shape[0] - drop_dupl_df.shape[0])
    )

    duplicated_gis = set(
        raw_df.groupby('gi_number', as_index=False) \
            .agg({'ass_id': 'count'}) \
            .query('ass_id > 1')['gi_number']
        )
    dupl_fpath = os.path.join(
        os.path.dirname(outfpath),
        'duplicated_GI_numbers.txt'
    )
    print('Writing duplicated GI numbers to file `{}`'.format(dupl_fpath))
    with open(dupl_fpath, 'w') as dupl_file:
        str_gis = map(str, duplicated_gis)
        dupl_file.write('\n'.join(str_gis) + '\n')
    # end with

    print('Removing duplicated lines in file `{}`'.format(outfpath))
    drop_dupl_df.to_csv(
        outfpath,
        sep='\t',
        index=False,
        header=True,
        encoding='utf-8',
        na_rep='NA',
        mode='w'
    )
else:
    print('No duplicated GI numbers found')
# end if

print('Completed!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
