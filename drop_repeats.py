#!/usr/bin/env python3

# Script removes gene sequences, which contain large repeats.
# Threshold (repeat length) for distinguishing "large" and "not large" repeats
#   is automatically defined.

# Input files:
# 1. `-i/--assm-acc-file` is output of script merge_assID2acc_and_remove_WGS.py.
#   It has 4 columns: ass_id, refseq_id, acc, title. `refseq_id` is GI number.
# 2. Fasta file with all extracted genes sequences (-f/--all-fasta-file).

# Output files:
# 1. Fasta file containing no sequences with NN (--out-fasta-file).
# 2. `--out-stats-file` is a TSV file containing per-replicon statisticsw for `--out-fasta-file`.
# 3. Fasta file containing sequences with NN (--NN-outfile).

import os
import sys
import argparse
from typing import Sequence, Dict

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from gene_seqs_2_stats import gene_seqs_2_stats


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--input-fasta-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-s',
    '--genes-stats-file',
    help='TSV file (with header) containing per-replicon SSU gene statistics',
    required=True
)

parser.add_argument(
    '-i',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '-r',
    '--repeats-file',
    help="""TSV file (with header) with
  coordinates and lengths of each repeat found by script find_repeats.py""",
    required=True
)

parser.add_argument(
    '-c',
    '--categories-file',
    help='per-gene categories file (output of script assign_genome_caregories.py)',
    required=True
)

# Output files

parser.add_argument(
    '--out-fasta',
    help='output fasta file containing sequences of genes without large repeats',
    required=True
)

parser.add_argument(
    '--out-stats',
    help='output per-replicon statistics file',
    required=True
)



args = parser.parse_args()


in_fasta_fpath = os.path.abspath(args.input_fasta_file)
in_stats_fpath = os.path.abspath(args.genes_stats_file)
assm_acc_fpath = os.path.abspath(args.assm_acc_file)
repeats_fpath = os.path.abspath(args.repeats_file)
cat_fpath = os.path.abspath(args.categories_file)
pure_genes_fpath = os.path.abspath(args.out_fasta)
pure_genes_stats_fpath = os.path.abspath(args.out_stats)


# Check existance of all input files
for fpath in (in_fasta_fpath, in_stats_fpath, assm_acc_fpath, repeats_fpath, cat_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directories if needed
for some_dir in map(os.path.dirname, [pure_genes_fpath, pure_genes_stats_fpath]):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end if


def select_gene_seqs(
    ass_id: str,
    seq_records: Sequence[SeqRecord],
    stats_df: pd.DataFrame) -> Dict[str, SeqRecord]:
    # Function selects SeqRecords from `seq_records` belonging to
    #   genome with Assembly ID `ass_id`.

    # Select ACCESSION.VERSIONs of current genome
    accs = set(stats_df[stats_df['ass_id'] == ass_id]['acc'])

    # Filter SeqRecords
    selected_seq_records = tuple(
        filter(
            lambda r: r.id.partition(':')[0] in accs,
            seq_records
        )
    )

    return {r.id: r for r in selected_seq_records}
# end def select_gene_seqs


def set_acc(row: pd.core.series.Series) -> pd.core.series.Series:
    # Function sets ACCESSION.VERSION to single row of a datraframe,
    #   derivind ACCESSION.VERSION from it's seqID.
    # This function is meant to be used by method pd.DataFrame.apply
    row['acc'] = row['seqID'].partition(':')[0]
    return row
# end def set_acc


def read_repeats_df(repeats_fpath: str) -> pd.DataFrame:
    # Function reads repeats dataframe and adds row indicating
    #   whether there is any conserved region in repeat.
    repeats_df = pd.read_csv(repeats_fpath, sep='\t')
    repeats_df['cons_reg_count'] = repeats_df['conserv_2'] \
                                 + repeats_df['conserv_3'] \
                                 + repeats_df['conserv_4'] \
                                 + repeats_df['conserv_5b'] \
                                 + repeats_df['conserv_6a'] \
                                 + repeats_df['conserv_7a'] \
                                 + repeats_df['conserv_8a'] \
                                 + repeats_df['conserv_9']

    return repeats_df
# end def read_repeats_df


# == Proceed ==

# Read file with Assembly IDs, ACCESSION.VERSIONs and titles
ass_acc_df = pd.read_csv(assm_acc_fpath, sep='\t')

# Read fle of per-replicon 16S statistics
in_stats_df = pd.read_csv(in_stats_fpath, sep='\t')

# Read categories file
cat_df = pd.read_csv(cat_fpath, sep='\t')

# Read sequence records
seq_records = tuple(SeqIO.parse(in_fasta_fpath, 'fasta'))

# Create set of input seqIDs
seqIDs = set( (r.id for r in seq_records) )

# Read repeats dataframe
repeats_df = read_repeats_df(repeats_fpath)

# Remove repeats, which do not contain conserved regions --
#   we will just keep them in result sequence file
repeats_df = repeats_df[repeats_df['cons_reg_count'] != 0]
repeats_df = repeats_df.query('seqID in @seqIDs')

# We need some additional info for finding repeat length threshold.
# Namely, we need Assembly IDs, categories and repeat lengths in the same dataframe.

repeats_df = repeats_df.merge(cat_df[['seqID', 'category', 'ass_id']], on='seqID', how='left') \
    .sort_values(by='rep_len', ascending=False) \
    .reset_index() \
    .drop_duplicates(subset=['seqID'], keep='first')


# == Find repeat length threshold ==

# We will start from the largest repeat and will drop all genes,
#   until we find either genome from non-3rd category or
#   genome which contains genes shorter than gene with repeat.

# Set initial threshold
repeat_len_threshold = repeats_df['rep_len'].max() + 1
if pd.isnull(repeat_len_threshold):
    repeat_len_threshold = float('inf')
# end if

print(f'\nmax repeat length = {repeat_len_threshold}')
print(repeats_df[['seqID', 'ass_id', 'gene_len', 'rep_len', 'category']].head(20))
print()

# Iterate over genes with repeats, starting from gene having the largest repeat
for i, row in repeats_df.iterrows():

    ass_id = row['ass_id']
    category = row['category']
    gene_len = row['gene_len']
    
    # Select gene sequence belonging to current genome
    selected_seq_records = select_gene_seqs(ass_id, seq_records, in_stats_df)

    # Check if shorter any gene exist in current genome
    genes_lengths = set([len(r.seq) for r in selected_seq_records.values()])
    shorter_exist = len(tuple(filter(lambda l: l < gene_len, genes_lengths))) > 0

    if shorter_exist or category == 3:
        repeat_len_threshold = row['rep_len'] # update threshold
    else:
        break
    # end if
# end for

# Subtract one from theshold in order not to miss "borderline" repeats --
#   those repeats, at which `repeat_len_threshold` wast updated last time.
# Below, we will use '>' operation to filter repeats by length.
repeat_len_threshold -= 1


# Print results of threshold detection
print(f'repeat_len_threshold = {repeat_len_threshold}')
print(repeats_df[repeats_df['rep_len'] <= repeat_len_threshold][['seqID', 'ass_id', 'gene_len', 'rep_len', 'category']].head(15))

# Read repets dataframe again, sincwe we modified it
repeats_df = read_repeats_df(repeats_fpath)
# Remove repeats, which do not contain conserved regions --
#   we will just keep them in result sequence file
repeats_df = repeats_df[repeats_df['cons_reg_count'] != 0]

# Remove seqIDs of sequences having large repeats, using detected threshold.
seqIDs_with_large_repeats = set(
    repeats_df[repeats_df['rep_len'] > repeat_len_threshold]['seqID']
)

# Remove sequences with long repeats
filtered_seq_records = tuple(
    filter(
        lambda r: r.id not in seqIDs_with_large_repeats,
        seq_records
    )
)

# Write filtered sequences to output file
print(f'Writing {len(filtered_seq_records)} seq records to file `{pure_genes_fpath}`...')
with open(pure_genes_fpath, 'wt') as fasta_outfile:
    for seq_record in filtered_seq_records:
        fasta_outfile.write(f'>{seq_record.description}\n{seq_record.seq}\n')
    # end for
# end with


print()

# Calculate per-replicon statistics
print('Calculating statistics')
gene_seqs_2_stats(pure_genes_fpath, assm_acc_fpath, pure_genes_stats_fpath)
print()

print('\nCompleted!')
print(pure_genes_fpath)
print(pure_genes_stats_fpath)
