#!/usr/bin/env python3

# Script maps seqIDs to TaxIDs using `--per-genome-taxID-file` file.

# Input files:
# 1. -i/--assm-acc-file -- is output of script merge_assID2acc_and_remove_WGS.py.
#   It has 4 columns: ass_id, gi_number, acc, title.
# 2. -f/--all-fasta-file -- fasta file of SSU gene sequences
# 3. --per-genome-taxID-file -- file mapping Assembly IDs to taxIDs

# Output files:
# 2. --per-gene-outfile -- output file mapping seqIDs to taxIDs

print(f'\n|=== STARTING SCRIPT `{__file__}` ===|\n')



import os
import sys
import argparse
from typing import Sequence


import pandas as pd
from Bio import SeqIO


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '-f',
    '--all-fasta-file',
    help='input fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '--per-genome-taxID-file',
    help='file mapping Assembly IDs to taxIDs',
    required=True
)

# Output files


parser.add_argument(
    '--per-gene-outfile',
    help='output file mapping genes seqIDs to taxIDs',
    required=True
)


args = parser.parse_args()


# For convenience
ass_acc_fpath = os.path.abspath(args.assm_acc_file)
fasta_seqs_fpath = os.path.abspath(args.all_fasta_file)
per_genome_taxID_fpath = os.path.abspath(args.per_genome_taxID_file)
per_gene_outfpath = os.path.abspath(args.per_gene_outfile)


# Check existance of all input files and dependencies
for fpath in (ass_acc_fpath, fasta_seqs_fpath, per_genome_taxID_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for


# Create output directories if needed
if not os.path.isdir(os.path.dirname(per_gene_outfpath)):
    try:
        os.makedirs(os.path.dirname(per_gene_outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(per_gene_outfpath)}`')
        sys.exit(1)
    # end try
# end if


def select_gene_seqIDs(ass_id: str,
    seq_records: Sequence[str],
    ass_acc_df: pd.DataFrame) -> Sequence[str]:

    # Get ACCESSION.VERSION's for current assembly
    accs = set(ass_acc_df[ass_acc_df['ass_id'] == ass_id]['acc'])

    # Filter genes from current genome
    selected_seq_records = tuple(
        filter(
            lambda r: r.id.partition(':')[0] in accs,
            seq_records
        )
    )

    # Make result dictionary and return it
    return tuple(map(lambda r: r.id, selected_seq_records))
# end def select_gene_seqs


ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')
per_genome_taxID_df = pd.read_csv(per_genome_taxID_fpath, sep='\t')

ass_acc_df = ass_acc_df.merge(
    per_genome_taxID_df[['ass_id', 'taxID']],
    on='ass_id',
    how='left'
)

seq_records = tuple(SeqIO.parse(fasta_seqs_fpath, 'fasta'))

ass_ids = set(ass_acc_df['ass_id'])


with open(per_gene_outfpath, 'wt') as per_gene_outfile:

    per_gene_outfile.write('ass_id\tseqID\ttaxID\n')

    for i, ass_id in enumerate(ass_ids):
        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        taxID = next(iter(
            set(ass_acc_df[ass_acc_df['ass_id'] == ass_id]['taxID'])
        ))

        curr_seqIDs = select_gene_seqIDs(ass_id, seq_records, ass_acc_df)

        for seqID in curr_seqIDs:
            per_gene_outfile.write(f'{ass_id}\t{seqID}\t{taxID}\n')
        # end for

    # end for
# end with


print('\nCompleted!')
print(per_gene_outfpath)
print(f'\n|=== EXITTING SCRIPT `{__file__}` ===|\n')
