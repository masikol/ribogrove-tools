#!/usr/bin/env python3

import os
import gzip

from Bio import SeqIO

accs_fpath = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/degenerate_in_16S/degenerate_in_16S_accs.txt'
gbk_dpath = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/genomes-dwnld/genomes-data/gbk'
outfpath = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/degenerate_in_16S/degenerate_in_16S_accs_and_NN.tsv'

accs = tuple(
    map(
        str.strip,
        open(accs_fpath).readlines()
    )
)

with open(outfpath, 'wt') as outfile:

    outfile.write('acc\tcontains_NN\n')

    for i, acc in enumerate(accs):
        print(f'\rDoing {i+1}/{len(accs)}: {acc}', end=' ')

        gbk_fpath = os.path.join(
            gbk_dpath,
            f'{acc}.gbk.gz'
        )

        with gzip.open(gbk_fpath, 'rt') as gbk_file:
            gbk_record = list(SeqIO.parse(gbk_file, 'gb'))[0]
        # end with

        contains_NN = 1 if 'NN' in str(gbk_record.seq).upper() else 0

        outfile.write(f'{acc}\t{contains_NN}\n')
    # end for
# end with

print('\nCompleted!')
print(outfpath)
