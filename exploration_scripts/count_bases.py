#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import sys
import argparse

from Bio import SeqIO

# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--input-fasta',
    help='input fasta file',
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
infpath = os.path.abspath(args.input_fasta)
outfpath = os.path.abspath(args.outfile)

# Check arguments
if not os.path.exists(infpath):
    print(f'Error: file `{fpath}` does not exist.')
    sys.exit(1)
# end if

if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create output directory `{os.path.dirname(outfpath)}`')
        print(str(err))
    # end try
# end if

print(infpath)
print()
# infpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/gene_seqs/bacteria_pure_gene_seqs_annotated.fasta'
# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/bases_count.tsv'


with open(infpath, 'rt') as infile, open(outfpath, 'w') as outfile:

    outfile.write(f'seqID\ta\tt\tg\tc\tr\ty\tw\ts\tk\tm\th\tv\tb\td\tn\tlen\n')

    for i, record in enumerate(SeqIO.parse(infile, 'fasta')):

        seq = str(record.seq).upper()

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

        outfile.write(f'{record.id}\t{a}\t{t}\t{g}\t{c}\t{r}\t{y}\t{w}\t{s}\t{k}\t{m}\t{h}\t{v}\t{b}\t{d}\t{n}\t{len(seq)}\n')
    # end for
# end with

print('\nCompleted!')
print(outfpath)
