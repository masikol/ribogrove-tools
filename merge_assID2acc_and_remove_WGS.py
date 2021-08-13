#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# Scripts merges TSV file, which is output of script assembly2refseq_id.py (`ass_gi_fpath`)
#   and TSV file, which is output of script gis_to_accs.py (`gi_acc_fpath`) on column `refseq_id`.
# Output (file `outfpath`) is a TSV file of 4 columns (ass_id, refseq_id, acc, title).

import numpy as np
import pandas as pd


ass_gi_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq.tsv'
gi_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs.tsv'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'


# Read input
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


def set_is_WGS(row):
    # Function sets `is_WGS` of a row to True if it's field `title` contains 
    #    "WHOLE GENOME SHOTGUN SEQUENCE". Otherwise sets it to False
    row['is_WGS'] = 'WHOLE GENOME SHOTGUN SEQUENCE' in row['title'].upper()
    return row
# end def set_is_WGS


# == Remove "WHOLE GENOME SHOTGUN SEQUENCE"s ==

print(f'len before rm WGS = {gi_acc_df.shape[0]}')

gi_acc_df['is_WGS'] = np.repeat(None, gi_acc_df.shape[0])
gi_acc_df = gi_acc_df.apply(set_is_WGS, axis=1)
gi_acc_df = gi_acc_df[gi_acc_df['is_WGS'] == False]

print(f'len after rm WGS = {gi_acc_df.shape[0]}')


# == Merge Assembly IDs with RefSeq Accession.Version ==

merged_df = ass_gi_df.merge(gi_acc_df, on='refseq_id', how='right')

print('MERGED: DATAFRAME')
print(merged_df.shape)
print(merged_df.head())

# Write output
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
