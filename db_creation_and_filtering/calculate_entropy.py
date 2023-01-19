#!/usr/bin/env python3

# The script calculates per-base intragenomic entropy from non-aberrant genes.
# The script aligns gene sequences with MUSCLE, and then caculates per-base entropy basing
#   on this multiple sequence alignment. The script calculates variability of the
#   target genes from the genomes of category 1 harbouring more than 1 target gene.

## Command line arguments
### Input files:
# 1. `-f / --fasta-seqs-file` -- a fasta file of gene sequences to be processed.
#   This is the file, which is the output of th script `make_final_seqs.py`. Mandatory.
# 2. `-s / --genes-stats-file` -- a TSV file of per-replicon genes statistics.
#   This is the file, which is the output of th script `make_final_seqs.py`. Mandatory.
# 3. `-c / --categories-file` -- a TSV file of categories info.
#   This is the file, which is the output of th script `assign_genome_categories.py`. Mandatory.

### Output files:
# 1. `-o / --outfile` -- a TSV file containing per-position intragenomic entropy. Mandatory.

### "Cached" files:
# 1. `--prev-per-base-entropy-file` -- a file of per-base entropy from the workdir
#   of the previous RiboGrove release (per_base_*_entropy.tsv.gz).
# 2. `--prev-assm-acc-file` -- a file `*_refseq_accs_merged.tsv` from the workdir
#   of the previous RiboGrove release.

### Dependencies:
# 1. `--muscle` -- a MUSCLE (https://www.drive5.com/muscle/) aligner executable.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')

import sys
import gzip
import math
import argparse
import operator
from array import array
import subprocess as sp
from io import StringIO
from functools import reduce
from typing import Sequence, Dict, List

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO



# == Parse arguments ==

parser = argparse.ArgumentParser()


# Input files

parser.add_argument(
    '-f',
    '--fasta-seqs-file',
    help='fasta file of SSU gene sequences without aberrant genes',
    required=True
)

parser.add_argument(
    '-s',
    '--genes-stats-file',
    help='TSV file (with header) containing per-replicons SSU gene statistics',
    required=True
)

parser.add_argument(
    '-c',
    '--categories-file',
    help='TSV file (with header) containing per-replicons SSU gene statistics',
    required=True
)

parser.add_argument(
    '--prev-per-base-entropy-file',
    help="""TSV file (with header: ass_id,pos,entropy)
    containing per-base intragenomic entropy. For example, file `per_base_demo_bacteria_entropy.tsv.gz`.
    The file may be gzipped.""",
    required=False
)

parser.add_argument(
    '--prev-assm-acc-file',
    help="""TSV file (with header) with
    Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs
    from previous RiboGrove release""",
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
    '--muscle',
    help='muscle executable',
    required=True
)


args = parser.parse_args()


# For convenience
fasta_seqs_fpath = os.path.abspath(args.fasta_seqs_file)
genes_stats_fpath = os.path.abspath(args.genes_stats_file)
categories_fpath = os.path.abspath(args.categories_file)
muscle_fpath = os.path.abspath(args.muscle)
outfpath = os.path.abspath(args.outfile)

cached_entropy = not args.prev_assm_acc_file         is None \
             and not args.prev_per_base_entropy_file is None
if cached_entropy:
    prev_perbase_entropy_fpath = os.path.abspath(args.prev_per_base_entropy_file)
    prev_assm_acc_fpath = os.path.abspath(args.prev_assm_acc_file)
else:
    cached_entropy = False
    prev_perbase_entropy_fpath = None
    prev_assm_acc_fpath = None
# end if


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, genes_stats_fpath, muscle_fpath, categories_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Check if executables are actually executable
if not os.access(muscle_fpath, os.X_OK):
    print(f'Error: file `{muscle_fpath}` is not executable!')
    sys.exit(1)
# end if

# Create output directory if needed
if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

# Check if previous ("cached") files is specified
if cached_entropy:
    for f in (prev_perbase_entropy_fpath, prev_assm_acc_fpath):
        if not os.path.exists(f):
            print(f'Error: file `{f}` does not exist!')
            sys.exit(1)
        # end if
    # end for
# end if


print(fasta_seqs_fpath)
print(genes_stats_fpath)
print(categories_fpath)
print(muscle_fpath)
print(prev_perbase_entropy_fpath)
print(prev_assm_acc_fpath)
print()


def select_gene_seqs(ass_id: str,
                     seq_records: Sequence[SeqRecord]) -> Sequence[SeqRecord]:
    selected_seq_records = filter(
        lambda r: get_ass_id_from_seq_record(r) == ass_id,
        seq_records
    )
    # Make result dictionary and return it
    return tuple(selected_seq_records)
# end def


def get_ass_id_from_seq_record(seq_record):
    return int(seq_record.id.partition(':')[0][2:])
# end def



def do_msa(seq_records: Sequence[SeqRecord], muscle_fpath: str) -> List[SeqRecord]:
    # Function does Multiple Sequence Alignment

    # Configure command
    cmd = f'{muscle_fpath} -quiet -diags'

    # Configure input fasta string for MSA
    fasta_str = '\n'.join(
        (f'>{r.id}\n{str(r.seq)}' for r in seq_records)
    ) + '\n'

    # Create pipe
    pipe = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    # Write inpit fasta record to stdin
    pipe.stdin.write(fasta_str.encode('ascii'))

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

    return msa_records
# end def do_msa


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

    return entropy_arr
# end def calc_entropy


def encode_accs(acc_list):
    return ''.join(sorted(acc_list))
# end def


def read_maybegzipped_df_tsv(fpath):
    if fpath.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    # end if

    with open_func(fpath, 'rt') as file:
        df = pd.read_csv(file, sep='\t')
    # end with

    return df
# end def


if cached_entropy:
    print('Reading cached files')

    prev_perbase_entropy_df = read_maybegzipped_df_tsv(prev_perbase_entropy_fpath)
    cached_ass_ids = set(prev_perbase_entropy_df['ass_id'])

    prev_assm_acc_df = pd.read_csv(prev_assm_acc_fpath, sep='\t')
    ass_ids = set(prev_assm_acc_df['ass_id'])
    acc_code_dict = dict()
    for ass_id in ass_ids:
        accs = tuple(prev_assm_acc_df[prev_assm_acc_df['ass_id'] == ass_id]['acc'])
        acc_code_dict[ass_id] = encode_accs(accs)
    # end for
    print('done')
else:
    cached_ass_ids = set()
# end if


# == Proceed ==

# Read statistics file
stats_df = pd.read_csv(genes_stats_fpath, sep='\t')

# Read categories file
categories_df = pd.read_csv(categories_fpath, sep='\t')

# Get Assembly IDs of 1 category
ass_ids = tuple(set(categories_df[categories_df['category'] == 1]['ass_id']))

# Read genes sequnces
seq_records = tuple(SeqIO.parse(fasta_seqs_fpath, 'fasta'))


per_base_entropy_fpath = os.path.join(
    os.path.dirname(outfpath),
    'per_base_' + os.path.basename(outfpath) + '.gz'
)

with gzip.open(per_base_entropy_fpath, 'wt') as per_base_entropy_outfile:

    per_base_entropy_outfile.write('ass_id\tpos\tentropy\n')

    # Iterate over assemblies
    for i, ass_id in enumerate(ass_ids):
        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        if ass_id in cached_ass_ids:
            acc_list = tuple(stats_df[stats_df['ass_id'] == ass_id]['acc'])
            acc_code = encode_accs(acc_list)
            if acc_code == acc_code_dict[ass_id]:
                cached_perbase_ass_df = prev_perbase_entropy_df[
                    prev_perbase_entropy_df['ass_id'] == ass_id
                ]
                cached_perbase_ass_df.to_csv(
                    per_base_entropy_outfile,
                    sep='\t',
                    header=False,
                    index=False,
                    na_rep='NA',
                    encoding='utf-8'
                )
                continue # cache hit
            # end if
        # end if

        # Select genes sequnences for currnet genome
        selected_seq_records = select_gene_seqs(ass_id, seq_records)

        # Perform MSA only if there are at least 2 sequences
        if len(selected_seq_records) > 1:
            # Perform MSA
            msa_records = do_msa(selected_seq_records, muscle_fpath)

            # Calculate entropy
            entropy_arr = calc_entropy(msa_records)

            # Write per-base entropy
            for i, entropy in enumerate(entropy_arr):
                per_base_entropy_outfile.write(f'{ass_id}\t{i}\t{entropy}\n')
            # end for
        # end if
    # end for
# end with

print()
print(f'The per-base entropy file is here: {per_base_entropy_fpath}')
print('Summarizing the calculated entropy...')


# Summarize entropy: calculate sum and mean for each genome

# Read per-base entropy file
with gzip.open(per_base_entropy_fpath, 'rt') as per_base_entropy_outfile:
    per_base_entropy_df = pd.read_csv(per_base_entropy_outfile, sep='\t')
# end with

# Calculate sum and mean entropy for each genome
# Calculate number of variable positions for each genome

count_var_positions = lambda entropy_arr: len(
    tuple(
        filter(
            lambda entropy: entropy > 1e-6,
            entropy_arr
        )
    )
)

summary_entropy_df = per_base_entropy_df.groupby('ass_id', as_index=False) \
    .agg({'entropy': ('sum', 'mean', count_var_positions)})

summary_entropy_df.columns = ['ass_id', 'sum_entropy', 'mean_entropy', 'num_var_cols']

del per_base_entropy_df

summary_entropy_df['num_var_cols'] = summary_entropy_df['num_var_cols'].map(int)


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
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
