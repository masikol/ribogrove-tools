#!/usr/bin/env python3

# The script calculates per-base intragenomic entropy from non-aberrant genes.
# The script aligns gene sequences with MAFFT, and then caculates per-base entropy basing
#   on this multiple sequence alignment. The script calculates variability of the
#   target genes from the genomes of category 1 harbouring more than 1 target gene.

## Command line arguments
### Input files:
# 1. `-f / --fasta-seqs-file` -- a fasta file of gene sequences to be processed.
#   This is the file, which is the output of th script `make_final_seqs.py`. Mandatory.
# 2. `-c / --categories-file` -- a TSV file of categories info.
#   This is the file, which is the output of th script `assign_genome_categories.py`. Mandatory.

### Parameters:
# 1. `-t / --threads` -- Number of CPU threads to use for MSA. Default: 1.
# 2. `--tmp-dir` -- A temporary directory for MAFFT-aligned sequences. Default: directory of the `--outfile`.

### Output files:
# 1. `-o / --outfile` -- a TSV file containing per-position intragenomic entropy. Mandatory.

### "Cached" files:
# 1. `--prev-per-base-entropy-file` -- a file of per-base entropy from the workdir
#   of the previous RiboGrove release (per_base_*_entropy.json.gz).

### Dependencies:
# 1. `--mafft` -- a MAFFT (https://mafft.cbrc.jp/alignment/software/) aligner executable.


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
    '-f',
    '--fasta-seqs-file',
    help='fasta file of SSU gene sequences without aberrant genes',
    required=True
)

parser.add_argument(
    '-c',
    '--categories-file',
    help='file mapping assembly accessions to genome categories',
    required=True
)

# Parameters

parser.add_argument(
    '--tmp-dir',
    help="""A temporary directory for MAFFT-aligned sequences.
Default: directory of the `--outfile`.""",
    required=False
)

parser.add_argument(
    '-t',
    '--threads',
    help='Number of CPU threads to use for MSA. Default: 1.',
    required=False,
    default=1
)

# "Cache"

parser.add_argument(
    '--prev-per-base-entropy-file',
    help="""a JSON file
    containing per-base intragenomic entropy. For example, file `per_base_demo_bacteria_entropy.json.gz`.
    The file may be gzipped.""",
    required=False
)

# Output files
parser.add_argument(
    '-o',
    '--outfile',
    help='output TSV file',
    required=True
)

# Dependencies
parser.add_argument(
    '--mafft',
    help='mafft executable',
    required=True
)

args = parser.parse_args()


# == Import them now ==
import sys
import gzip
import math
import json
import operator
from array import array
import subprocess as sp
from io import StringIO
from functools import reduce
from typing import Sequence, Dict, List

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.ribogrove_seqID import parse_asm_acc


# For convenience
fasta_seqs_fpath = os.path.abspath(args.fasta_seqs_file)
categories_fpath = os.path.abspath(args.categories_file)
mafft_fpath = os.path.abspath(args.mafft)
outfpath = os.path.abspath(args.outfile)
if args.tmp_dir is None:
    tmp_dirpath = os.path.dirname(outfpath)
else:
    tmp_dirpath = os.path.abspath(args.tmp_dir)
# end if

try:
    threads = int(args.threads)
    if threads < 1:
        raise ValueError
    # end if
except ValueError:
    print('Error: threads must be an integer > 0. Your value: `{}`'.format(threads))
# end try

cache_mode = not args.prev_per_base_entropy_file is None
if cache_mode:
    prev_perbase_entropy_fpath = os.path.abspath(args.prev_per_base_entropy_file)
else:
    prev_perbase_entropy_fpath = None
# end if


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, mafft_fpath, categories_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist')
        sys.exit(1)
    # end if
# enb for

# Check if executables are actually executable
if not os.access(mafft_fpath, os.X_OK):
    print(f'Error: file `{mafft_fpath}` is not executable')
    sys.exit(1)
# end if

# Create output and tmp directories if needed
for d in (tmp_dirpath, os.path.dirname(outfpath)):
    if not os.path.isdir(d):
        try:
            os.makedirs(d)
        except OSError as err:
            print(f'Error: cannot create directory `{d}`')
            sys.exit(1)
        # end try
    # end if
# enf for

# Check if previous ("cached") files is specified
if cache_mode:
    if not os.path.exists(prev_perbase_entropy_fpath):
        print(f'Error: file `{prev_perbase_entropy_fpath}` does not exist')
        sys.exit(1)
    # end if
# end if


print(fasta_seqs_fpath)
print(categories_fpath)
print(mafft_fpath)
print(threads)
print(tmp_dirpath)
print(prev_perbase_entropy_fpath)
print()

def select_gene_seqs(asm_acc: str,
                     seq_records: Sequence[SeqRecord]) -> Sequence[SeqRecord]:
    selected_seq_records = filter(
        lambda r: get_asm_acc_from_seq_record(r) == asm_acc,
        seq_records
    )
    # Make result dictionary and return it
    return tuple(selected_seq_records)
# end def


def get_asm_acc_from_seq_record(seq_record):
    return parse_asm_acc(seq_record.id)
# end def



def do_msa(seq_records: Sequence[SeqRecord],
           mafft_fpath: str,
           threads : int,
           tmp_dirpath : str) -> List[SeqRecord]:
    # Function does Multiple Sequence Alignment

    # TODO: remove hardcoded threads
    tmp_fpath = os.path.join(tmp_dirpath, 'tmp.fasta')

    # Configure command
    cmd = ' '.join(
        [
            mafft_fpath,
            '--auto',
            '--thread 2',
            tmp_fpath
        ]
    )

    # Configure input fasta string for MSA
    with open(tmp_fpath, 'wt') as outfile:
        for seq_record in seq_records:
            outfile.write('>{}\n{}\n'.format(seq_record.id, seq_record.seq))
    # end with

    # Create pipe
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)

    # Run command
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error while doing msa!')
        print('seqIDs:')
        print(' '.join( [r.id for r in seq_records] ))
        print(stdout_stderr[1].decode('utf-8'))
    else:
        # Get StringIO of result alignment
        msa_io = StringIO(stdout_stderr[0].decode('utf-8'))
    # end if

    # Parse alignment
    msa_records = list(SeqIO.parse(msa_io, 'fasta'))
    msa_io.close()

    os.unlink(tmp_fpath)

    return msa_records
# end def


def get_aln_column(i: int, seqs: Sequence[str]) -> Sequence[str]:
    # Function extracts single column from multiple sequence alignment
    return tuple(map(lambda s: s[i], seqs))
# end


def calc_entropy(msa_records: Sequence[SeqRecord]) -> Sequence[float]:
    # Function calculates intragenomic per-base entropy.

    n_seqs = len(msa_records)

    # Extract aligned sequences as type `str`
    seqs = tuple(map(lambda x: str(x.seq), msa_records))
    seq_length = len(seqs[0]) # they are aligned and have the same length

    # Init entropy array
    entropy_arr = array('d', np.repeat(np.nan, seq_length))
    aln_pos = 0

    # Iterate over positions in alignment
    for i in range(seq_length):

        # Extract alignment column of index i
        aln_column = get_aln_column(i, seqs)

        # Count frequences aof all bases at this position
        freqs_arr = tuple(
            map(
                lambda base: aln_column.count(base) / n_seqs,
                set(aln_column)
            )
        )

        # Calculate entropy
        # abs instead of minus in order not to allow "-0.0" values
        entropy_arr[aln_pos] = abs(
            reduce(
                operator.add,
                (freq * math.log(freq, 2) for freq in freqs_arr)
            )
        )
        aln_pos += 1
    # end for

    return tuple(entropy_arr)
# end def


def encode_accs(acc_list):
    return ''.join(sorted(acc_list))
# end def


def set_sum_mean_num_var_cols(row):
    global perbase_entropy_dict
    asm_acc = row['asm_acc']
    entropy_arr = perbase_entropy_dict[asm_acc]
    row['sum_entropy'] = np.sum(entropy_arr)
    row['mean_entropy'] = np.mean(entropy_arr)
    row['num_var_cols'] = count_var_positions(entropy_arr)
    return row
# end def

def count_var_positions(entropy_arr):
    return len(
        tuple(
            filter(
                lambda entropy: entropy > 1e-6,
                entropy_arr
            )
        )
    )
# end def


# == Proceed ==

# Read categories file
categories_df = pd.read_csv(categories_fpath, sep='\t')

# Get Assembly IDs of 1 category
asm_accs = frozenset(
    categories_df[categories_df['category'] == 1]['asm_acc']
)

# Read genes sequnces
seq_records = tuple(SeqIO.parse(fasta_seqs_fpath, 'fasta'))


per_base_entropy_fpath = os.path.join(
    os.path.dirname(outfpath),
    'per_base_entropy.json.gz'
)

if cache_mode:
    print('Creating auxiliary data structures...')

    with gzip.open(prev_perbase_entropy_fpath, 'rt') as input_handle:
        prev_perbase_entropy_dict = json.load(input_handle)
    # end with
    prev_asm_accs = frozenset(prev_perbase_entropy_dict.keys())

    cached_asm_accs = asm_accs & prev_asm_accs
    del prev_asm_accs
    print('done')
    print(
        '{:,}/{:,} genomes are cached' \
            .format(len(cached_asm_accs), len(asm_accs))
    )
else:
    prev_perbase_entropy_dict = dict()
    cached_asm_accs = set()
# end if



perbase_entropy_dict = {asm_acc: None for asm_acc in asm_accs}

# Iterate over assemblies
for i, asm_acc in enumerate(asm_accs):
    print(f'\rDoing {i+1}/{len(asm_accs)}: {asm_acc}', end=' '*10)
    if asm_acc in cached_asm_accs:
        perbase_entropy_dict[asm_acc] = prev_perbase_entropy_dict[asm_acc]
        continue # cache hit
    # end if

    # Select genes sequnences for currnet genome
    selected_seq_records = select_gene_seqs(asm_acc, seq_records)

    # Perform MSA only if there are at least 2 sequences
    if len(selected_seq_records) > 1:
        # Perform MSA
        msa_records = do_msa(
            selected_seq_records,
            mafft_fpath,
            threads,
            tmp_dirpath
        )
        # Calculate entropy
        perbase_entropy_dict[asm_acc] = calc_entropy(msa_records)
    else:
        del perbase_entropy_dict[asm_acc]
    # end if
# end for

# Save per-base entropy
with gzip.open(per_base_entropy_fpath, 'wt') as output_handle:
    json.dump(perbase_entropy_dict, output_handle)
# end with

print()
print(f'The per-base entropy file is here: `{per_base_entropy_fpath}`')
print('Summarizing the calculated entropy...')


# Summarize entropy: calculate sum and mean for each genome

# Calculate sum and mean entropy for each genome
# Calculate number of variable positions for each genome


if len(perbase_entropy_dict) != 0:
    summary_entropy_df = pd.DataFrame(
        {
            'asm_acc': list(perbase_entropy_dict.keys()),
            'sum_entropy': np.repeat(None, len(perbase_entropy_dict)),
            'mean_entropy': np.repeat(None, len(perbase_entropy_dict)),
            'num_var_cols': np.repeat(None, len(perbase_entropy_dict)),
        }
    )
    summary_entropy_df = summary_entropy_df.apply(set_sum_mean_num_var_cols, axis=1)
else:
    summary_entropy_df = pd.DataFrame(
        {
            'asm_acc': [],
            'sum_entropy': [],
            'mean_entropy': [],
            'num_var_cols': [],
        }
    )
# end if


# Overwrite the per-base file: it is barely informative
summary_entropy_df.to_csv(
    outfpath,
    sep='\t',
    mode='w',
    index=False,
    header=True
)


print('\nCompleted!')
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
