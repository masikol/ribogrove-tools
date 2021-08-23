#!/usr/bin/env python3

# Script calculates per-base intragenomic entropy from non-aberrant genes.

# Input files:
# 1. Fasta file of genes sequences (-f/--fasta-seqs-file).
# 2. TSV file if per-replicon genes statistics (-s/--genes-stats-file).

# Output files:
# 1. TSV file containing per-position intragenomic entropy (-o/--outfile).

# Dependencies:
# 1. MUSCLE aligner (--muscle).


import os
import sys
import math
import argparse
import operator
from array import array
import subprocess as sp
from io import StringIO
from functools import reduce
from typing import Sequence, List

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
muscle_fpath = os.path.abspath(args.muscle)
outfpath = os.path.abspath(args.outfile)


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, genes_stats_fpath, muscle_fpath):
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


def select_gene_seqs(ass_id: str,
    seq_records: Sequence[str],
    stats_df: pd.DataFrame) -> Sequence[SeqRecord]:

    # Get ACCESSION.VERSION's for current assembly
    accs = set(stats_df[stats_df['ass_id'] == ass_id]['acc'])

    # Filter genes from current genome
    selected_seq_records = tuple(
        filter(
            lambda r: r.id.partition(':')[0] in accs,
            seq_records
        )
    )

    # Make result dictionary and return it
    return selected_seq_records
# end def select_gene_seqs



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


# == Proceed ==

# Read statistics file
stats_df = pd.read_csv(genes_stats_fpath, sep='\t')

# Get Assembly IDs
ass_ids = tuple(set(stats_df['ass_id']))

# Read genes sequnces
seq_records = tuple(SeqIO.parse(fasta_seqs_fpath, 'fasta'))


with open(outfpath, 'wt') as entropy_outfile:

    entropy_outfile.write('ass_id\tpos\tentropy\n')

    # Iterate over assemblies
    for i, ass_id in enumerate(ass_ids):
        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        # Select genes sequnences for currnet genome
        selected_seq_records = select_gene_seqs(ass_id, seq_records, stats_df)

        # Perform MSA only if there are at least 2 sequences
        if len(selected_seq_records) > 1:

            # Perform MSA
            msa_records = do_msa(selected_seq_records, muscle_fpath)

            # Calculate entropy
            entropy_arr = calc_entropy(msa_records)

            # Write per-base entropy
            for i, entropy in enumerate(entropy_arr):
                entropy_outfile.write(f'{ass_id}\t{i}\t{entropy}\n')
            # end for
        # end if
    # end for
# end with


print('\nCompleted!')
print(outfpath)
