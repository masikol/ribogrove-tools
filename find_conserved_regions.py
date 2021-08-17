#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils


input_fasta_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'
stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_stats_no_NN.tsv'
conserved_regions_fpath = '/mnt/1.5_drive_0/16S_scrubbling/consensus_seqs/conserved_regions_NR.fasta'

outfpath = '/mnt/1.5_drive_0/16S_scrubbling/consensus_seqs/conserved_regions_locations.tsv'


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


conserved_seq_records = tuple(SeqIO.parse(conserved_regions_fpath, 'fasta'))


stats_df = pd.read_csv(stats_fpath, sep='\t')
ass_ids = set(stats_df['ass_id'])



seq_records = tuple(SeqIO.parse(input_fasta_fpath, 'fasta'))


with open(outfpath, 'wt') as outfile:

    outfile.write('ass_id\tseqID\tcons_region\tstart\n')

    for i, ass_id in enumerate(ass_ids):

        print(f'\r Doing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        selected_seq_records = select_gene_seqs(ass_id, seq_records, stats_df)

        for seqID, seq_record in selected_seq_records.items():
            for conserved_seq_record in conserved_seq_records:
                search_result = SeqUtils.nt_search(str(seq_record.seq), str(conserved_seq_record.seq))
                if len(search_result) > 1:
                    for start in search_result[1:]:
                        outfile.write(f'{ass_id}\t{seqID}\t{conserved_seq_record.id}\t{start+1}\n')
                    # end for
                # end if
            # end for

        # end for


    # end for
# end with


print('\nCompleted!')
print(outfpath)
