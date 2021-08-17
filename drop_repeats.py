#!/usr/bin/env python3

import os

import numpy as np
import pandas as pd
from Bio import SeqIO

from gene_seqs_2_stats import gene_seqs_2_stats


# no_aberrant_genes_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/non_aberrant_genes.fasta'
no_aberrant_genes_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'
# non_aberrant_genes_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/non_aberrant_genes_stats.tsv'
non_aberrant_genes_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_stats_no_NN.tsv'

ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
repeats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/repeats_no_NN.tsv'
cat_fpath = '/mnt/1.5_drive_0/16S_scrubbling/categories/bacteria_genome_categories.tsv'

pure_genes_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/test_pure_genes_seqs.fasta'
pure_genes_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/test_pure_genes_stats.tsv'

# pure_genes_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/pure_genes_seqs.fasta'
# pure_genes_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/pure_genes_stats.tsv'


def select_gene_seqs(ass_id, seq_records, stats_df):

    accs = set(stats_df[stats_df['ass_id'] == ass_id]['acc'])

    selected_seq_records = tuple(
        filter(
            lambda r: r.id.partition(':')[0] in accs,
            seq_records
        )
    )

    return {r.id: r for r in selected_seq_records}
# end def select_gene_seqs


def set_acc(row):
    row['acc'] = row['seqID'].partition(':')[0]
    return row
# end def set_acc


def read_repeats_df(repeats_fpath):
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


os.system(f'seqkit stats -a {no_aberrant_genes_fpath}')


ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')

non_aberrant_genes_stats_df = pd.read_csv(non_aberrant_genes_stats_fpath, sep='\t')

cat_df = pd.read_csv(cat_fpath, sep='\t')


repeats_df = read_repeats_df(repeats_fpath)
repeats_df = repeats_df[repeats_df['cons_reg_count'] != 0]

# We need some additional info for finding repeat length threshold.
# Namely, we need assembly IDs, categories and repeat lengths in the same dataframe.

repeats_df['acc'] = np.repeat(None, repeats_df.shape[0])
repeats_df = repeats_df.apply(set_acc, axis=1)

repeats_df = ass_acc_df[['ass_id', 'acc']].merge(repeats_df, on='acc', how='right')
repeats_df = repeats_df.groupby('ass_id').agg({'rep_len': 'max'}) \
    .reset_index() \
    .merge(cat_df[['ass_id', 'category']], on='ass_id', how='left') \
    .sort_values(by='rep_len', ascending=False) \
    .reset_index()


seq_records = tuple(SeqIO.parse(no_aberrant_genes_fpath, 'fasta'))


# == Find repeat length thresgold ==

repeat_len_threshold = repeats_df['rep_len'].max() + 1

print(f'max repeat length = {repeat_len_threshold}')

for i, row in repeats_df.iterrows():

    ass_id = row['ass_id']
    category = row['category']
    
    selected_seq_records = select_gene_seqs(ass_id, seq_records, non_aberrant_genes_stats_df)

    genes_are_identical = len(set([str(r.seq) for r in selected_seq_records.values()])) < 2
    # genes_are_identical = len(set([len(r.seq) for r in selected_seq_records.values()])) < 2

    if not genes_are_identical or category == 3:
        repeat_len_threshold = row['rep_len']
    else:
        break
    # end if
# end for

repeat_len_threshold -= 1


print(f'repeat_len_threshold = {repeat_len_threshold}')


# Read repets dataframe again
repeats_df = read_repeats_df(repeats_fpath)
repeats_df = repeats_df[repeats_df['cons_reg_count'] != 0]

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

# Write filtered sequences to the "pure" file
print(f'Writing {len(filtered_seq_records)} seq records to file `{pure_genes_fpath}`...')
with open(pure_genes_fpath, 'wt') as fasta_outfile:
    for seq_record in filtered_seq_records:
        fasta_outfile.write(f'>{seq_record.description}\n{seq_record.seq}\n')
    # end for
# end with


print()
os.system(f'seqkit stats -a {pure_genes_fpath}')

print('Calculating statistics')
gene_seqs_2_stats(pure_genes_fpath, ass_acc_fpath, pure_genes_stats_fpath)
print()

print('\nCompleted!')
print(pure_genes_fpath)
print(pure_genes_stats_fpath)
