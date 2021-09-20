#!/usr/bin/env python3

import os
import sys
import gzip
import statistics as sts

import pandas as pd
from Bio import SeqIO

in_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
gbk_dpath = '/mnt/1.5_drive_0/16S_scrubbling/genomes-data/gbk'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/all_genes_notes.tsv'

# in_acc_fpath = '/home/deynonih/cager/new_16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
# gbk_dpath = '/home/deynonih/Desktop/gbk'
# outfpath = '/home/deynonih/cager/new_16S_scrubbling/all_genes_notes.tsv'


# note:
# 16S ribosomal RNA rRNA prediction is too short
# note:
# possible 16S ribosomal RNA but does not have goodblast hits on one or both of the ends


acc_df = pd.read_csv(
    in_acc_fpath,
    sep='\t'
)

needed_types = {'RRNA', 'MISC_FEATURE'}

n_accs = acc_df.shape[0]

with open(outfpath, 'wt') as outfile:

    outfile.write('ass_id\tacc\tfeature_type\tfeature_start\tnote\n')

    for i, row in acc_df.iterrows():
        ass_id = row['ass_id']
        acc = row['acc']
        print(f'\rDoing {i+1}/{n_accs}: {acc}', end=' '*10)

        gbk_fpath = os.path.join(
            gbk_dpath,
            f'{acc}.gbk.gz'
        )

        with gzip.open(gbk_fpath, 'rt') as gbfile:
            gbrecord = tuple(SeqIO.parse(gbfile, 'gb'))[0]
        # end with

        for feature in gbrecord.features:
            if feature.type.upper() in needed_types:
                qualifiers = feature.qualifiers
                if 'note' in qualifiers.keys():
                    for note in qualifiers['note']:
                        feature_start = feature.location.start
                        outfile.write(f'{ass_id}\t{acc}\t{feature.type}\t{feature_start}\t{note}\n')
                    # end if
                # end if
            # end if
        # end for

        # print('-' * 20)
        # input()

        # fasta_outfile.write(f'{acc}\t{seqtech}\n')
    # end for
# end with

print('\nCompleted!')
print(outfpath)
