#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import time

from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'

# Setup
assm_id_fpath = '/mnt/1.5_drive_0/16S_scrubbling/assembly_UIDs.txt'
# assm_id_fpath = sys.argv[1]
if not os.path.exists(assm_id_fpath):
    print('Error: file `{}` does not exist'.format(assm_id_fpath))
    sys.exit(0)
# end if

outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq.tsv'
# assm_id_fpath = sys.argv[2]


ass_ids = tuple(
    map(
        str.strip,
        open(assm_id_fpath).readlines()
    )
)

# done_ass_ids = set(
#     map(
#         lambda l: l.split('\t')[0],
#         open(outfpath).readlines()
#     )
# )

# with open(outfpath, 'at') as outfile:
with open(outfpath, 'wt') as outfile:

    outfile.write(f'ass_id\trefseq_id\n')

    for i, ass_id in enumerate(ass_ids):

        # if ass_id in done_ass_ids:
        #     continue
        # # end if

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' ')

        n_errors = 0

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

        try:
            refseq_ids = [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']]
        except IndexError:
            pass
        else:
            for refseq_id in refseq_ids:
                outfile.write(f'{ass_id}\t{refseq_id}\n')
            # end for
        # end try

        time.sleep(0.4)
    # end for
# end with

print('\nCompleted!')
print(outfpath)
