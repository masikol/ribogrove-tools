#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# The script takes output of script assembly2refseq_id.py (RefSeq GI numbers) as input (file `infpath`)
#   and translates them to "accession.version"s and titles.
# Output (file `outfpath`) is a TSV file of 3 columns (refseq_id, acc, title).


import os
import sys
import time

import pandas as pd

from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'

infpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq.tsv'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs.tsv'


# Read input
gi_df = pd.read_csv(
    infpath,
    sep='\t'
)


chunk_size = 50
n_done_ids = 0

print('\r0/{}'.format(gi_df.shape[0]), end=' ')


# == Proceed ==
with open(outfpath, 'wt') as outfile:

    # Write header
    outfile.write(f'refseq_id\tacc\ttitle\n')

    # Iterate over chunks of RefSeq IDs and get their "accession.version"s and titles
    for i in range(0, gi_df.shape[0], chunk_size):

        curr_gis = tuple(
            map(
                str,
                tuple(gi_df.iloc[i : i + chunk_size,]['refseq_id'])
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
            except:
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

print('\nCompleted!')
print(outfpath)
