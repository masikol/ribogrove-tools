# -*- encoding: utf-8 -*-

import os
import sys
from io import StringIO
# import subprocess as sp
import statistics as sts

import pandas as pd
from Bio import SeqIO


def _select_seqs(acc, seq_records):

    selected_seq_records = tuple(
        filter(
            lambda r: acc in r.id,
            seq_records
        )
    )

    return selected_seq_records
# end def _select_seqs


def _get_len(record):
    return len(record.seq)
# end def _get_len


def gene_seqs_2_stats(seqs_fasta_fpath, ass_acc_fpath, stats_outfpath):

    ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')

    seq_records = tuple(SeqIO.parse(seqs_fasta_fpath, 'fasta'))

    with open(stats_outfpath, 'wt') as stats_outfile:

        stats_outfile.write('ass_id\trefseq_id\tacc\ttitle\tnum_genes\tmin_len\tmax_len\tmean_len\tmedian_len\n')

        for i, row in ass_acc_df.iterrows():

            ass_id = row['ass_id']
            refseq_id = row['refseq_id']
            acc = row['acc']
            title = row['title']

            print(f'\rDoing {i+1}/{ass_acc_df.shape[0]}: {acc}', end=' '*10)

            selected_seq_records = _select_seqs(acc, seq_records)

            num_genes = len(selected_seq_records)

            if num_genes != 0:
                gene_lengths = tuple(map(_get_len, selected_seq_records))
                min_len = min(gene_lengths)
                max_len = max(gene_lengths)
                mean_len = sts.mean(gene_lengths)
                median_len = sts.median(gene_lengths)
            else:
                min_len = 'NA'
                max_len = 'NA'
                mean_len = 'NA'
                median_len = 'NA'
            # end if

            stats_outfile.write(f'{ass_id}\t{refseq_id}\t{acc}\t{title}\t{num_genes}\t')
            stats_outfile.write(f'{min_len}\t{max_len}\t{mean_len}\t{median_len}\n')

        # end for
    # end with
# end def gene_seqs_2_stats
