#!/usr/bin/env python3

import os
import glob

import pandas as pd


acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
outdir = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/genomes-dwnld/genomes-data/gbk'

accs = tuple(
    pd.read_csv(
        acc_fpath,
        sep='\t'
    )['acc']
)

num_not_exist = 0
num_prev_version = 0

for acc_version in accs:
    acc = acc_version.partition('.')[0]
    outfpath = os.path.join(outdir, f'{acc_version}.gbk.gz')
    if not os.path.exists(outfpath):
        num_not_exist += 1
        if len(glob.glob(os.path.join(outdir, f'{acc}*'))) != 0:
            num_prev_version += 1
        # end if
    # end if
# end for

print(f'{num_not_exist}/{len(accs)}')
print(f'{num_prev_version}/{num_not_exist}')
