#!/usr/bin/env python3

import os
import sys
import time
import http
import gzip
from functools import partial

import pandas as pd
from Bio import Entrez


Entrez.email = "maximdeynonih@gmail.com"

acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
outdir = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/genomes-dwnld/genomes-data/gbk'
log_fpath = '/mnt/1.5_drive_0/16S_scrubbling/genomes-download.0.log'


def write_to_file(text, log_fpath):
    if text != '':
        with open(log_fpath, 'a') as logfile:
            logfile.write(f'{text}\n')
        # end with
    # end if
# end def write_to_file

write_to_log_file = partial(write_to_file, log_fpath=log_fpath)


with open(log_fpath, 'w') as _:
    pass
# end with

n_done = 0

accs = tuple(
    pd.read_csv(
        acc_fpath,
        sep='\t'
    )['acc']
)


n_accs = len(accs)
max_errors = 3
downloaded = 0

for i, acc in enumerate(accs):

    print(f'\rDoing {i+1}/{n_accs}: {acc}', end=' '*10)

    outfpath = os.path.join(outdir, f'{acc}.gbk.gz')

    if os.path.exists(outfpath):
        write_to_log_file(f'{acc} - ok (already exists)')
        continue
    # end if

    downloaded += 1

    errors = 0
    while errors < max_errors:

        try:
            gbk_handle = Entrez.efetch(
                db='nucleotide',
                id=acc,
                # rettype='gb',
                rettype='gbwithparts',
                retmode='text',
                linkname='refseq'
            )

            with gzip.open(outfpath, 'wt') as outfile:
                outfile.write(gbk_handle.read())
            # end with

            gbk_handle.close()

        except IOError as err:
            err_msg = f'Error: cannot fetch GenBank file for `{acc}`: {err}'
            print(f'\n{err_msg}')
            write_to_log_file(err_msg)
        except (RuntimeError, http.client.IncompleteRead):
            errors += 1
        else:
            break
        # end try

        time.sleep(0.4)
    # end while

    if errors >= max_errors:
        err_msg = f'Error: 3 times RuntimeError occured for `{acc}`'
        print(f'\n{err_msg}')
        write_to_log_file(err_msg)
    else:
        write_to_log_file(f'{acc} - ok')
    # end if

    print(f'\r{i+1}/{n_accs} done', end=' '*10)

# end for

print('\n\nCompleted!')
print(f'downloaded = {downloaded}')
