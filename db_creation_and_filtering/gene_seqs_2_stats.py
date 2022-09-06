# -*- encoding: utf-8 -*-

# Module contains funciton `gene_seqs_2_stats`.
# The function creates per-replicon statistics from given fasta file containing gene sequences.
# See it's descrition below (Ctrl+F "def gene_seqs_2_stats").


import statistics as sts

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def _get_len(record: SeqRecord):
    # Function returns length of given SeqRecord `record`
    return len(record.seq)
# end def _get_len


def gene_seqs_2_stats(seqs_fasta_fpath: str, ass_acc_fpath: str, stats_outfpath: str) -> None:
    # Function creates per-replicon statistics from given fasta file containing gene sequences.
    # :param seqs_fasta_fpath: path to input fasta file;
    # :param ass_acc_fpath: path to file containing columns (ass_id, gi_number, acc, title);
    #   (it is output of script merge_assID2acc_and_remove_WGS.py)
    # :param stats_outfpath: path to output per-replicon statistics file;


    # Read Assembly-ACCESSION.VERSION file
    ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')

    # Read seq_reords from input fasta file
    seq_records = SeqIO.parse(seqs_fasta_fpath, 'fasta')

    print('Creating auxiliary data structures...')
    seq_record_dict = dict()
    for seq_record in seq_records:
        acc = seq_record.id.split(':')[1]
        try:
            seq_record_dict[acc].append(seq_record)
        except KeyError:
            seq_record_dict[acc] = [seq_record]
        # end try
    # end for
    print('Done')

    # == Proceed ==

    print(f'0/{ass_acc_df.shape[0]} replicons processed', end=' '*10)

    with open(stats_outfpath, 'wt') as stats_outfile:

        # Write header
        stats_outfile.write('ass_id\tgi_number\tacc\ttitle\tnum_genes\tmin_len\tmax_len\tmean_len\tmedian_len\n')

        status_bar_step = 500
        next_update = 500

        # Iterate over row of Assembly-ACCESSION.VERSION dataframe
        for i, row in ass_acc_df.iterrows():

            ass_id = row['ass_id']
            gi_number = row['gi_number']
            acc = row['acc']
            title = row['title']

            # Select SeqRecords for current replicon
            # selected_seq_records = _select_seqs(acc, seq_records)
            try:
                selected_seq_records = seq_record_dict[acc]
            except KeyError:
                num_genes = 0 # no genes for this replicon
            else:
                # Tally selected genes
                num_genes = len(selected_seq_records)
                del seq_record_dict[acc]
            # end if


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
            stats_outfile.write(f'{ass_id}\t{gi_number}\t{acc}\t{title}\t{num_genes}\t')
            stats_outfile.write(f'{min_len}\t{max_len}\t{mean_len}\t{median_len}\n')

            if i + 1 > next_update:
                print(f'\r{i+1}/{ass_acc_df.shape[0]} replicons processed', end=' '*10)
                next_update += status_bar_step
            # end if

        # end for
        print(f'\r{ass_acc_df.shape[0]}/{ass_acc_df.shape[0]} replicons processed')
    # end with
# end def gene_seqs_2_stats
