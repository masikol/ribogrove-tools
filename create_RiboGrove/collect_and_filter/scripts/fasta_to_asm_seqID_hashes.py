#!/usr/bin/env python3

import os
import sys
import argparse

from src.rg_tools_time import get_time


sys.stderr.write(
    '\n|=== {} STARTING SCRIPT `{}` ===|\n\n'.format(
        get_time(), os.path.basename(__file__)
    )
)


parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input-fasta', required=True)
parser.add_argument('-o', '--output-tsv', required=True)
args = parser.parse_args()


# == Import them now ==

import gzip

from Bio import SeqIO
from src.util import make_asm_seqID_hash


in_fasta_fpath = os.path.abspath(args.input_fasta)
out_tsv_fpath = os.path.abspath(args.output_tsv)

if not os.path.exists(in_fasta_fpath):
    sys.stderr.write('Error: input file `{}` does not exist!\n'.format(in_fasta_fpath))
    sys.exit(1)
# end if

out_dir = os.path.dirname(out_tsv_fpath)
if out_dir and not os.path.isdir(out_dir):
    try:
        os.makedirs(out_dir)
    except OSError as err:
        sys.stderr.write('Error: cannot create directory `{}`\n'.format(out_dir))
        sys.exit(1)
    # end try
# end if


asm_to_seqIDs: dict[str, set[str]] = {}

open_func = gzip.open if in_fasta_fpath.endswith('.gz') else open
with open_func(in_fasta_fpath, 'rt') as infile:
    for record in SeqIO.parse(infile, 'fasta'):
        seqID = record.id
        asm_acc = seqID.partition(':')[0]
        asm_to_seqIDs.setdefault(asm_acc, set()).add(seqID)
    # end for
# end with

with open(out_tsv_fpath, 'wt') as outfile:
    outfile.write('asm_acc\tasm_seqID_hash\n')
    for asm_acc in asm_to_seqIDs:
        seqIDs_sorted = asm_to_seqIDs[asm_acc]
        hash_val = make_asm_seqID_hash(seqIDs_sorted)
        outfile.write('{}\t{}\n'.format(asm_acc, hash_val))
    # end for
# end with

sys.stderr.write(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n\n'.format(
        get_time(), os.path.basename(__file__)
    )
)

sys.exit(0)
