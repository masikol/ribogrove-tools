#!/usr/bin/env python3

import os
import sys
import subprocess as sp
from io import StringIO

import pandas as pd
from Bio import SeqIO

from gene_seqs_2_stats import gene_seqs_2_stats


bioawk = '/home/cager/Misc_soft/bioawk-1.0/bioawk'

raw_genes_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'
aberrant_seqIDs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/aberrant_seqIDs.txt'
ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
categories_fpath = '/mnt/1.5_drive_0/16S_scrubbling/categories/bacteria_16S_genes_categories.tsv'

non_aberrant_genes_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/non_aberrant_genes.fasta'
non_aberrant_genes_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/non_aberrant_genes_stats.tsv'


def remove_aberrant_genes(raw_genes_seqs_fpath, aberrant_seqIDs_fpath, non_aberrant_genes_seqs_fpath):

    cmd = f'cat {raw_genes_seqs_fpath} | seqkit grep -vf {aberrant_seqIDs_fpath} > {non_aberrant_genes_seqs_fpath}'
    print(cmd)
    os.system(cmd)
# end def remove_aberrant_genes


def remove_3cat_sort_tail(non_aberrant_genes_seqs_fpath, categories_fpath):

    categories_df = pd.read_csv(categories_fpath, sep='\t')

    get_length_cmd = f'cat {non_aberrant_genes_seqs_fpath} | {bioawk} -c fastx ' \
        + '\'{print $name "\\t" length($seq)}\''

    pipe = sp.Popen(get_length_cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Eror while extracting gene lengths')
        print(stdout_stderr[1].decode('utf-8'))
        print('cmd:')
        print(f'  {get_length_cmd}')
        sys.exit(1)
    else:
        gene_lengths_str = 'seqID\tlen\n' + stdout_stderr[0].decode('utf-8')
    # end if

    gene_lengths_io = StringIO(gene_lengths_str)
    lengths_df = pd.read_csv(gene_lengths_io, sep='\t') \
        .sort_values(by='len', ascending=True).reset_index()
    gene_lengths_io.close()

    lengths_df = lengths_df.merge(categories_df, on='seqID', how='left')

    first_non_cat3_index = 0
    for i, row in lengths_df.iterrows():
        if row['category'] != 3:
            first_non_cat3_index = i
            break
        # end if
    # end for

    print(lengths_df.shape)
    print(lengths_df.head())

    print(f'Cutting {first_non_cat3_index} shortest genes')
    lengths_df = lengths_df.iloc[first_non_cat3_index:]
    print(f'{lengths_df.shape[0]} genes left')
    print(f'Now, min length is {lengths_df["len"].min()}')

    print(lengths_df.shape)
    print(lengths_df.head())

    seqIDs_to_keep = set(lengths_df['seqID'])

    seq_records = tuple(SeqIO.parse(non_aberrant_genes_seqs_fpath, 'fasta'))

    no_3cat_tail_seq_records = tuple(
        filter(
            lambda r: r.id in seqIDs_to_keep,
            seq_records
        )
    )

    print(f'Writing {len(seq_records)} seq records to `{non_aberrant_genes_seqs_fpath}`...')
    with open(non_aberrant_genes_seqs_fpath, 'wt') as non_aberrant_genes_seqs_outfile:
        for record in no_3cat_tail_seq_records:
            non_aberrant_genes_seqs_outfile.write(f'>{record.description}\n{str(record.seq)}\n')
        # end for
    # end with
# end def remove_3cat_sort_tail


os.system(f'seqkit stats -a {raw_genes_seqs_fpath}')

remove_aberrant_genes(raw_genes_seqs_fpath, aberrant_seqIDs_fpath, non_aberrant_genes_seqs_fpath)

os.system(f'seqkit stats -a {non_aberrant_genes_seqs_fpath}')

remove_3cat_sort_tail(non_aberrant_genes_seqs_fpath, categories_fpath)
# remove_3cat_sort_tail(raw_genes_seqs_fpath, categories_fpath)

os.system(f'seqkit stats -a {non_aberrant_genes_seqs_fpath}')


print('Calculating statistics')
gene_seqs_2_stats(non_aberrant_genes_seqs_fpath, ass_acc_fpath, non_aberrant_genes_stats_fpath)
print()

print('\nCompleted!')
print(non_aberrant_genes_stats_fpath)
print(non_aberrant_genes_seqs_fpath)
