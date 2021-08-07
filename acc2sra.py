#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import time

import pandas as pd
from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'

# Setup
acc_fpath = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/degenerate_in_16S/degenerate_in_16S_accs_and_NN.tsv'
# acc_fpath = sys.argv[1]
if not os.path.exists(acc_fpath):
    print('Error: file `{}` does not exist'.format(acc_fpath))
    sys.exit(0)
# end if

outfpath = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/degenerate_in_16S/degenerate_in_16S_accs_2_sra.tsv'
# acc_fpath = sys.argv[2]


acc_df = pd.read_csv(
    acc_fpath,
    sep='\t'
)
acc_df = acc_df[acc_df['contains_NN'] == 0].reset_index()
# acc_df = acc_df.iloc[0:5,]

# done_ass_ids = set(
#     map(
#         lambda l: l.split('\t')[0],
#         open(outfpath).readlines()
#     )
# )

# with open(outfpath, 'at') as outfile:
with open(outfpath, 'wt') as outfile:

    outfile.write(f'acc\tbiosample_ids\tsra_ids\n')

    for i, row in acc_df.iterrows():
    # for i, ass_id in enumerate(ass_ids):

        biosample_ids = None
        sra_ids = None

        acc = row['acc']

        # if ass_id in done_ass_ids:
        #     continue
        # # end if

        print(f'\rDoing {i+1}/{acc_df.shape[0]}: {acc}', end=' ')

        n_errors = 0

        while n_errors < 3:
            try:
                handle = Entrez.elink(
                    dbfrom='nuccore',
                    db='biosample',
                    id=acc
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
            biosample_ids = [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']]
        except IndexError:
            biosample_ids = list()
        # else:
        #     for refseq_id in refseq_ids:
        #         outfile.write(f'{ass_id}\t{refseq_id}\n')
        #     # end for
        # end try

        time.sleep(0.5)

        if len(biosample_ids) == 0:
            outfile.write(f'{acc}\tNA\tNA\n')
            continue
        # end if


        n_errors = 0

        while n_errors < 3:
            try:
                handle = Entrez.elink(
                    dbfrom='biosample',
                    db='sra',
                    id=','.join(biosample_ids)
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
            sra_ids = [link['Id'] for link in record[0]['LinkSetDb'][0]['Link']]
        except IndexError:
            sra_ids = list()
        # end try

        print(sra_ids)

        time.sleep(0.5)

        outfile.write('{}\t{}\t{}\n'.format(
            acc,
            ";".join(biosample_ids) if len(biosample_ids) != 0 else 'NA',
            ";".join(sra_ids) if len(sra_ids) != 0 else 'NA'
        ))
    # end for
# end with

print('\nCompleted!')
print(outfpath)
