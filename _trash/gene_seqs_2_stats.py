# -*- encoding: utf-8 -*-

# Module contains funciton `gene_seqs_2_stats`.
# The function creates per-replicon statistics from given fasta file containing gene sequences.
# See it's descrition below (Ctrl+F "def gene_seqs_2_stats").


import statistics as sts

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.ribogrove_seqID import parse_seq_acc


def gene_seqs_2_stats(seqs_fasta_fpath: str, replicon_map_fpath: str, stats_outfpath: str) -> None:
    # Function creates per-replicon statistics from given fasta file containing gene sequences.
    # :param seqs_fasta_fpath: path to input fasta file;
    # :param replicon_map_fpath: TODO: describe param
    # :param stats_outfpath: path to output per-replicon statistics file;


    # Read Assembly-ACCESSION.VERSION file
    replicon_map_df = pd.read_csv(replicon_map_fpath, sep='\t')

    # Read seq_reords from input fasta file
    seq_records = SeqIO.parse(seqs_fasta_fpath, 'fasta')

    print('Creating auxiliary data structures...')
    seq_record_dict = dict()
    for seq_record in seq_records:
        seq_acc = parse_seq_acc(seq_record.id)
        try:
            seq_record_dict[seq_acc].append(seq_record)
        except KeyError:
            seq_record_dict[seq_acc] = [seq_record]
        # end try
    # end for
    print('Done')

    # == Proceed ==

    print(f'0/{replicon_map_df.shape[0]} replicons processed', end=' '*10)

    with open(stats_outfpath, 'wt') as stats_outfile:

        # Write header
        stats_outfile.write('asm_acc\tseq_acc\tnum_genes\n')

        status_bar_step = 500
        next_update = 500

        # Iterate over row of Assembly-ACCESSION.VERSION dataframe
        for i, row in replicon_map_df.iterrows():

            asm_acc = row['asm_acc']
            seq_acc = row['seq_acc']

            # Select SeqRecords for current replicon
            try:
                selected_seq_records = seq_record_dict[seq_acc]
            except KeyError:
                num_genes = 0 # no genes for this replicon
            else:
                # Tally selected genes
                num_genes = len(selected_seq_records)
                del seq_record_dict[seq_acc]
            # end if

            # Calculate statistics if at least one gene of current replicon
            #   is in `seq_records`
            gene_lengths = tuple(map(_get_len, selected_seq_records))

            # Write result line to output file
            stats_outfile.write(f'{asm_acc}\t{seq_acc}\t{num_genes}\n')

            if i + 1 > next_update:
                print(f'\r{i+1}/{replicon_map_df.shape[0]} replicons processed', end=' '*10)
                next_update += status_bar_step
            # end if

        # end for
        print(f'\r{replicon_map_df.shape[0]}/{replicon_map_df.shape[0]} replicons processed')
    # end with
# end def


def _get_len(record: SeqRecord):
    # Function returns length of given SeqRecord `record`
    return len(record.seq)
# end def
