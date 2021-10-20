#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# The Script takes IDs (IDs of database NCBI Assembly) from file `-i/--assm-id-file` (one per line)
#   and translates them to RefSeq GI numbers using elink utility (https://www.ncbi.nlm.nih.gov/books/NBK25497/).
# The script writes result to TSV file `-o/--outfile` of 2 columns (ass_ID, refseq_id)


import os
import sys
import time
import argparse


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
    outfile.write('ass_id\trefseq_id\n')

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
            except IOError as err:
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
            refseq_ids = [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']]
        except IndexError:
            pass
        else:
            # Write output
            for refseq_id in refseq_ids:
                outfile.write(f'{ass_id}\t{refseq_id}\n')
            # end for
        # end try

        # Wait a bit: we don't want NCBI to ban us :)
        time.sleep(0.4)
    # end for
# end with

print('\nCompleted!')
print(outfpath)
