#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import gzip

import pandas as pd
from Bio import SeqIO


gbk_dir_path = '/mnt/1.5_drive_0/16S_scrubbling/genomes-data/gbk'

# ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/bacteria_ass_refseq_accs_merged.tsv'
# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/bases_count_whole_genomes.tsv'

ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/archaea_refseq_accs_merged.tsv'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/bases_count_whole_genomes.tsv'


DEGEN_BASES = ('R', 'Y', 'W', 'S', 'K', 'M', 'H', 'V', 'B', 'D', 'N')


ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')

ass_ids = set(ass_acc_df['ass_id'])


def retrieve_seq(acc):
    gbk_fpath = os.path.join(gbk_dir_path, f'{acc}.gbk.gz')
    with gzip.open(gbk_fpath, 'rt') as gbk_file:
        gbk_record = tuple(SeqIO.parse(gbk_file, 'gb'))[0]
    # end with

    return str(gbk_record.seq).upper()
# end def retrieve_seq


def remove_degen_bases(seq):
    global DEGEN_BASES

    for b in DEGEN_BASES:
        seq = seq.replace(b, '')
    # end for

    return seq
# end def remove_degen_bases


with open(outfpath, 'w') as outfile:

    outfile.write(f'ass_id\ta\tt\tg\tc\ts\tlen\tlen_no_degen\n')

    for i, ass_id in enumerate(ass_ids):

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        accs = tuple(ass_acc_df[ass_acc_df['ass_id'] == ass_id]['acc'])

        seq = ''.join(
            map(
                retrieve_seq,
                accs
            )
        )

        init_len = len(seq)

        seq = remove_degen_bases(seq)

        s = seq.count('S')

        len_no_degen = len(seq)

        a = seq.count('A')
        t = seq.count('T')
        g = seq.count('G')
        c = seq.count('C')

        outfile.write(f'{ass_id}\t{a}\t{t}\t{g}\t{c}\t{s}\t{init_len}\t{len_no_degen}\n')
    # end for
# end with

print('\nCompleted!')
print(outfpath)
