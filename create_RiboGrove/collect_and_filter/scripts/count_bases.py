#!/usr/bin/env python3

# The script counts bases in the sequences in a fasta file:
#   bases `A`, `T`, `G`, `C`, `R`, `Y`, `W`, `S`, `K`, `M`, `H`, `V`, `B`, `D`, `N`.

## Command line arguments
### Input files:
# 1. `-i / --input-fasta` -- an input fasta file.

### Output files:
# 1. `-o / --outfile` -- an output TSV file.


import os
from src.rg_tools_time import get_time

print(
    '\n|=== {} STARTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)


# == Parse arguments ==
import argparse

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


# == Import them now ==
import sys

from Bio import SeqIO


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
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
