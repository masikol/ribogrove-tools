#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# Script takes IDs (IDs of database NCBI Assembly) from file `assm_id_fpath` (one per line)
#   and translates them to RefSeq GI numbers using elink utility (https://www.ncbi.nlm.nih.gov/books/NBK25497/).
# The script writes result to TSV file `outfpath` of 2 columns (ass_ID, refseq_id)


import os
import sys
import time

from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'


assm_id_fpath = '/mnt/1.5_drive_0/16S_scrubbling/assembly_UIDs.txt'

outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq.tsv'


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
    outfile.write(f'ass_id\trefseq_id\n')

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
