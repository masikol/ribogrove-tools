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
    print(f'Error: file `{infpath}` does not exist.')
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

SEP = '\t'

with open(infpath, 'rt') as infile, open(outfpath, 'w') as outfile:

    outfile.write('{}\n'.format(SEP.join([
        'seqID',
        'a', 't', 'g', 'c',
        'r', 'y', 'w', 's', 'k', 'm',
        'h', 'v', 'b', 'd',
        'n',
        'len',
    ])))

    for i, record in enumerate(SeqIO.parse(infile, 'fasta')):

        seq = str(record.seq).upper()

        a = str(seq.count('A'))
        t = str(seq.count('T'))
        g = str(seq.count('G'))
        c = str(seq.count('C'))
        r = str(seq.count('R'))
        y = str(seq.count('Y'))
        w = str(seq.count('W'))
        s = str(seq.count('S'))
        k = str(seq.count('K'))
        m = str(seq.count('M'))
        h = str(seq.count('H'))
        v = str(seq.count('V'))
        b = str(seq.count('B'))
        d = str(seq.count('D'))
        n = str(seq.count('N'))

        outfile.write('{}\n'.format(SEP.join([
            record.id,
            a, t, g, c,
            r, y, w, s, k, m,
            h, v, b, d,
            n,
            str(len(seq))
        ])))
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

sys.exit(0)
