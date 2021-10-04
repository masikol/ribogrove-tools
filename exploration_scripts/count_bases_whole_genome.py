#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import sys
import gzip
import argparse

import pandas as pd
from Bio import SeqIO

# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-a',
    '--ass-acc-file',
    help='input TSV file of Assembly IDs mapped t o ACCESSION.VERSIONs',
    required=True
)

parser.add_argument(
    '-g',
    '--gbk-dir',
    help='directory of input `.gbk.gz` files',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='output TSV file',
    required=True
)

args = parser.parse_args()

# For convenience
ass_acc_fpath = os.path.abspath(args.ass_acc_file)
gbk_dir_path = os.path.abspath(args.gbk_dir)
outfpath = os.path.abspath(args.outfile)


# gbk_dir_path = '/mnt/1.5_drive_0/16S_scrubbling/genomes-data/gbk'

# ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/bacteria_refseq_accs_merged.tsv'
# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/bases_count_whole_genomes.tsv'

# ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/archaea_refseq_accs_merged.tsv'
# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/bases_count_whole_genomes.tsv'


if not os.path.exists(ass_acc_fpath):
    print(f'Error: input file `{ass_acc_fpath}` does not exist.')
    sys.exit(1)
# end if

if not os.path.isdir(gbk_dir_path):
    print(f'Error: input directory `{gbk_dir_path}` does not exist.')
    sys.exit(1)
# end if

if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Cannot create output directory `{os.path.dirname(outfpath)}`')
        print(str(err))
        sys.exit(1)
    # end try
# end if


# DEGEN_BASES = ('R', 'Y', 'W', 'S', 'K', 'M', 'H', 'V', 'B', 'D', 'N')


ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')

ass_ids = set(ass_acc_df['ass_id'])


def retrieve_seq(acc):
    gbk_fpath = os.path.join(gbk_dir_path, f'{acc}.gbk.gz')
    with gzip.open(gbk_fpath, 'rt') as gbk_file:
        gbk_record = tuple(SeqIO.parse(gbk_file, 'gb'))[0]
    # end with

    return str(gbk_record.seq).upper()
# end def retrieve_seq


# def remove_degen_bases(seq):
#     global DEGEN_BASES

#     for b in DEGEN_BASES:
#         seq = seq.replace(b, '')
#     # end for

#     return seq
# # end def remove_degen_bases


with open(outfpath, 'w') as outfile:

    outfile.write(f'ass_id\ta\tt\tg\tc\tr\ty\tw\ts\tk\tm\th\tv\tb\td\tn\tlen\tlen_no_degen\n')

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

        a = seq.count('A')
        t = seq.count('T')
        g = seq.count('G')
        c = seq.count('C')
        r = seq.count('R')
        y = seq.count('Y')
        w = seq.count('W')
        s = seq.count('S')
        k = seq.count('K')
        m = seq.count('M')
        h = seq.count('H')
        v = seq.count('V')
        b = seq.count('B')
        d = seq.count('D')
        n = seq.count('N')

        # seq = remove_degen_bases(seq)
        # len_no_degen = len(seq)

        # outfile.write(f'{ass_id}\t{a}\t{t}\t{g}\t{c}\t{r}\t{y}\t{w}\t{s}\t{k}\t{m}\t{h}\t{v}\t{b}\t{d}\t{n}\t{init_len}\t{len_no_degen}\n')
        outfile.write(f'{ass_id}\t{a}\t{t}\t{g}\t{c}\t{r}\t{y}\t{w}\t{s}\t{k}\t{m}\t{h}\t{v}\t{b}\t{d}\t{n}\t{init_len}\n')
    # end for
# end with

print('\nCompleted!')
print(outfpath)
