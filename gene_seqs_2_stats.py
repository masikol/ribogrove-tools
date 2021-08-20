# -*- encoding: utf-8 -*-

# Module contains funciton `gene_seqs_2_stats`.
# The function creates per-replicon statistics from given fasta file containing gene sequences.
# See it's descrition below (Ctrl+F "def gene_seqs_2_stats").


import os
import sys
from io import StringIO
import statistics as sts
from typing import List, Sequence

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _select_seqs(acc: str, seq_records: Sequence[SeqRecord]) -> Sequence[SeqRecord]:
    # Function selects SeqRecords of given ACCESSION.VERSION (`acc`).

    selected_seq_records = tuple(
        filter(
            lambda r: acc in r.id,
            seq_records
        )
    )

    return selected_seq_records
# end def _select_seqs


def _get_len(record: SeqRecord):
    # Function returns length of given SeqRecord `record`
    return len(record.seq)
# end def _get_len


def gene_seqs_2_stats(seqs_fasta_fpath: str, ass_acc_fpath: str, stats_outfpath: str) -> None:
    # Function creates per-replicon statistics from given fasta file containing gene sequences.
    # :param seqs_fasta_fpath: path to input fasta file;
    # :param ass_acc_fpath: path to file containing columns (ass_id, refseq_id, acc, title);
    #   (it is output of script merge_assID2acc_and_remove_WGS.py)
    # :param stats_outfpath: path to output per-replicon statistics file;


    # Read Assembly-ACCESSION.VERSION file
    ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')

    # Read seq_reords from input fasta file
    seq_records = tuple(SeqIO.parse(seqs_fasta_fpath, 'fasta'))

    # == Proceed ==
    with open(stats_outfpath, 'wt') as stats_outfile:

        # Write header
        stats_outfile.write('ass_id\trefseq_id\tacc\ttitle\tnum_genes\tmin_len\tmax_len\tmean_len\tmedian_len\n')

        # Iterate over row of Assembly-ACCESSION.VERSION dataframe
        for i, row in ass_acc_df.iterrows():

            ass_id = row['ass_id']
            refseq_id = row['refseq_id']
            acc = row['acc']
            title = row['title']

            print(f'\rDoing {i+1}/{ass_acc_df.shape[0]}: {acc}', end=' '*10)

            # Select SeqRecords for current replicon
            selected_seq_records = _select_seqs(acc, seq_records)

            # Tally selected genes
            num_genes = len(selected_seq_records)

            # Calculate statistics if at least one gene of current replicon
            #   is in `seq_records`
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

            # Write result line to output file
            stats_outfile.write(f'{ass_id}\t{refseq_id}\t{acc}\t{title}\t{num_genes}\t')
            stats_outfile.write(f'{min_len}\t{max_len}\t{mean_len}\t{median_len}\n')
        # end for
    # end with
# end def gene_seqs_2_stats
