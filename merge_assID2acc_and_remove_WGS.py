#!/usr/bin/env python3

import numpy as np
import pandas as pd

ass_gi_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq.tsv'
gi_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs.tsv'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'

ass_gi_df = pd.read_csv(
    ass_gi_fpath,
    sep='\t'
)

gi_acc_df = pd.read_csv(
    gi_acc_fpath,
    sep='\t'
)

print(ass_gi_df.shape)
print(ass_gi_df.head())

print(gi_acc_df.shape)
print(gi_acc_df.head())


# Remove 'whole genome shotgun sequence'

def set_is_WGS(row):
    row['is_WGS'] = 'WHOLE GENOME SHOTGUN SEQUENCE' in row['title'].upper()
    return row
# end def set_is_WGS

print(f'len before rm WGS = {gi_acc_df.shape[0]}')

gi_acc_df['is_WGS'] = np.repeat(None, gi_acc_df.shape[0])
gi_acc_df = gi_acc_df.apply(set_is_WGS, axis=1)
gi_acc_df = gi_acc_df[gi_acc_df['is_WGS'] == False]

print(f'len after rm WGS = {gi_acc_df.shape[0]}')

# Merge Assembly IDs with RefSeq Accession.Version

merged_df = ass_gi_df.merge(gi_acc_df, on='refseq_id', how='right')

print('MERGED: DATAFRAME')
print(merged_df.shape)
print(merged_df.head())

merged_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    na_rep='NA',
    encoding='utf-8',
    columns=['ass_id', 'refseq_id', 'acc', 'title',]
)

print('\nCompleted!')
print(outfpath)
