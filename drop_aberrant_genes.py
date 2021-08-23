#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess as sp
from io import StringIO

import pandas as pd
from Bio import SeqIO

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
    '-a',
    '--aberrant-genes-file',
    help='file of seqIDs of aberrant genes, one per line',
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
    '-c',
    '--categories-file',
    help="""TSV file (with header) mapping seqIDs to categories""",
    required=True
)


# Output files

parser.add_argument(
    '--output-fasta',
    help='output fasta file which will contain non-aberrant genes\' sequences',
    required=True
)

parser.add_argument(
    '--output-stats',
    help='output TSV file of per-replicon SSU genes statistics',
    required=True
)


# Dependencies

parser.add_argument(
    '--seqkit',
    help='seqkit executable',
    required=True
)


args = parser.parse_args()


# For convenience
input_genes_seqs_fpath = os.path.abspath(args.input_fasta_file)
aberrant_seqIDs_fpath = os.path.abspath(args.aberrant_genes_file)
ass_acc_fpath = os.path.abspath(args.assm_acc_file)
categories_fpath = os.path.abspath(args.categories_file)
non_aberrant_genes_seqs_fpath = os.path.abspath(args.output_fasta)
non_aberrant_genes_stats_fpath = os.path.abspath(args.output_stats)
seqkit_fpath = os.path.abspath(args.seqkit)


def remove_aberrant_genes(
    input_genes_seqs_fpath: str,
    aberrant_seqIDs_fpath: str,
    non_aberrant_genes_seqs_fpath: str) -> None:
    # Function removes aberrant genes using seqkit grep

    # Configure command:
    # 1) read from `input_genes_seqs_fpath`
    # 2) remove sequences having seqIDs stored in file `aberrant_seqIDs_fpath`, one per line
    # 3) write remaining sequences to `non_aberrant_genes_seqs_fpath`
    cmd = f'cat {input_genes_seqs_fpath} | seqkit grep -vf {aberrant_seqIDs_fpath} > {non_aberrant_genes_seqs_fpath}'

    print(cmd)

    # Run command
    os.system(cmd)
# end def remove_aberrant_genes


def remove_3cat_sort_tail(non_aberrant_genes_seqs_fpath: str, categories_fpath: str) -> None:
    # Function removes so called "3-category short tail" from non-aberrant genes.
    # These are those genes, which are short and obviously aberrant and originate
    #   from 3-category genome.
    # We will sort genes by length ascendingly. Then, starting from the shortest gene,
    #   we will remove all genes until we reach first non-3-category gene.

    # Read per-gene categories dataframe
    categories_df = pd.read_csv(categories_fpath, sep='\t')

    # Read input genes sequences
    non_aberrant_seq_records = tuple(SeqIO.parse(non_aberrant_genes_seqs_fpath, 'fasta'))

    # Create length dataframe (seqID, len)
    lengths_df = pd.DataFrame(
        {
            'seqID': tuple(
                map(
                    lambda r: r.id, non_aberrant_seq_records
                )
            ),
            'len':  tuple(
                map(
                    lambda r: len(r.seq), non_aberrant_seq_records
                )
            ),
        }
    )

    del non_aberrant_seq_records

    # Sort genes by length in ascending order
    lengths_df = lengths_df.sort_values(by='len', ascending=True).reset_index()

    # Add categories to seqIDs and lengths
    lengths_df = lengths_df.merge(categories_df, on='seqID', how='left')

    # Find first non-3-category gene.
    first_non_cat3_index = 0
    for i, row in lengths_df.iterrows():
        if row['category'] != 3:
            first_non_cat3_index = i
            break
        # end if
    # end for

    print(f'Cutting {first_non_cat3_index} shortest genes')

    print(lengths_df.shape)
    print(lengths_df.head())

    # Remove "3-category short tail"
    lengths_df = lengths_df.iloc[first_non_cat3_index:]

    print(f'{lengths_df.shape[0]} genes left')
    print(f'Now, min length is {lengths_df["len"].min()}')

    print(lengths_df.shape)
    print(lengths_df.head())


    # Now we will actually filter genes
    seqIDs_to_keep = set(lengths_df['seqID'])

    # Read fil
    seq_records = tuple(SeqIO.parse(non_aberrant_genes_seqs_fpath, 'fasta'))

    # Filter genes: exclude those in "3-category short tail"
    no_3cat_tail_seq_records = tuple(
        filter(
            lambda r: r.id in seqIDs_to_keep,
            seq_records
        )
    )

    # Write filtered sequences to the save file
    print(f'Writing {len(seq_records)} seq records to `{non_aberrant_genes_seqs_fpath}`...')
    with open(non_aberrant_genes_seqs_fpath, 'wt') as non_aberrant_genes_seqs_outfile:
        for record in no_3cat_tail_seq_records:
            non_aberrant_genes_seqs_outfile.write(f'>{record.description}\n{str(record.seq)}\n')
        # end for
    # end with
# end def remove_3cat_sort_tail


# == Proceed ==

# Remove aberrant genes discovered by script `find_aberrant_genes.py`
remove_aberrant_genes(input_genes_seqs_fpath, aberrant_seqIDs_fpath, non_aberrant_genes_seqs_fpath)

# Remove "3-category short tail"
remove_3cat_sort_tail(non_aberrant_genes_seqs_fpath, categories_fpath)

print(non_aberrant_genes_seqs_fpath)


# Calculate per-replicon statistics of filtered genes

print('Calculating statistics')
gene_seqs_2_stats(non_aberrant_genes_seqs_fpath, ass_acc_fpath, non_aberrant_genes_stats_fpath)
print()

print('\nCompleted!')
print(non_aberrant_genes_stats_fpath)
