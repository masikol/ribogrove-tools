#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from io import StringIO

import pandas as pd
from Bio import SeqIO


# in_fasta_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'
in_fasta_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/pure_genes_seqs.fasta'

tax_fpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/per_gene_taxonomy.tsv'
cat_fpath = '/mnt/1.5_drive_0/16S_scrubbling/categories/bacteria_16S_genes_categories.tsv'

# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/TEST_gene_seqs_no_NN_annotated.fasta'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/pure_genes_seqs_annotated.fasta'

tax_sep = ';'


tax_df = pd.read_csv(
    tax_fpath,
    sep='\t',
    dtype={
        'seqID': str,
        'ass_id': pd.Int32Dtype(),
        'accs': str,
        'taxID': pd.Int32Dtype(),
        'tax_name': str,
        'genus': str,
        'family': str,
        'order': str,
        'class': str,
        'phylum': str,
        'superkingdom': str,
    }
)

cat_df = pd.read_csv(
    cat_fpath,
    sep='\t',
    dtype={
        'ass_id': pd.Int32Dtype(),
        'seqID': str,
        'category': pd.Int8Dtype(),
    }
)

n_seqs = len(tuple(SeqIO.parse(in_fasta_fpath, 'fasta')))

tax_df.index = tax_df['seqID']


with open(outfpath, 'wt') as outfile:

    seq_records = SeqIO.parse(in_fasta_fpath, 'fasta')

    for i, seq_record in enumerate(seq_records):

        print(f'\r Doing {i+1}/{n_seqs}: {seq_record.id}', end=' '*10)

        curr_tax_record = tax_df.loc[seq_record.id, ]

        taxonomy = tax_sep.join(
            (
                'Bateria',
                'NA' if pd.isnull(curr_tax_record['phylum']) else curr_tax_record['phylum'],
                'NA' if pd.isnull(curr_tax_record['class']) else curr_tax_record['class'],
                'NA' if pd.isnull(curr_tax_record['order']) else curr_tax_record['order'],
                'NA' if pd.isnull(curr_tax_record['family']) else curr_tax_record['family'],
                'NA' if pd.isnull(curr_tax_record['genus']) else curr_tax_record['genus'],
            )
        )

        if not pd.isnull(curr_tax_record['tax_name']):
            tax_name = curr_tax_record['tax_name'].replace(' ', '_')
        else:
            tax_name = 'no_taxonomy_name'
        # end if

        category = cat_df[cat_df['seqID'] == seq_record.id]['category'].values[0]

        seq_record.description = f'{seq_record.id} {tax_name} {tax_sep}{taxonomy}{tax_sep} category:{category}'
        outfile.write(f'>{seq_record.description}\n{seq_record.seq}\n')
    # end for
# end with

tmp_fasta_fpath = './tmp.fasta'
os.system(f'cat {outfpath} | seqkit seq -u > {tmp_fasta_fpath}')
os.system(f'cat {tmp_fasta_fpath} > {outfpath}')
os.unlink(tmp_fasta_fpath)

os.system(f'seqkit stats -a {outfpath}')

print('\nCompleted!')
print(outfpath)
