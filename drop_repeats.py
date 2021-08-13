#!/usr/bin/env python3

import os

from Bio import SeqIO

from gene_seqs_2_stats import gene_seqs_2_stats


max_repeat_len = 23

ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
repeats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/repeats_no_NN.tsv'
no_aberrant_genes_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/non_aberrant_genes.fasta'
non_aberrant_genes_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/non_aberrant_genes_stats.tsv'
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


os.system(f'seqkit stats -a {no_aberrant_genes_fpath}')


ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')
repeats_df = pd.read_csv(repeats_fpath, sep='\t')
cat_df = pd.read_csv(cat_fpath, sep='\t')
non_aberrant_genes_stats_df = pd.read_csv(non_aberrant_genes_stats_fpath, sep='\t')


cat_dict = {
    row['ass_id']: row['category'] \
        for i, row in cat_df.iterrows()
}

gene_count_dict = {
    row['ass_id']: row['num_genes'] \
        for i, row in non_aberrant_genes_stats_df.groupby('ass_id').agg({'num_genes': 'sum'}) \
            .reset_index().iterrows()
}

seqIDs_with_large_repeats = set(
    repeats_df[repeats_df['rep_len'] > max_repeat_len]['seqID']
)


seq_records = SeqIO.parse(no_aberrant_genes_fpath, 'fasta')


# ass_ids = set(ass_acc_df['ass_id'])
ass_ids = {4513671, 4513681, 7294851}

with open(pure_genes_fpath, 'wt') as fasta_outfile:
    for i, ass_id in enumerate(ass_ids):

        print(f'Doing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        selected_seq_records = select_gene_seqs(ass_id, seq_records, non_aberrant_genes_stats_df)

        if len(selected_seq_records) == 0:
            continue
        # end if

        seqIDs = set(map(lambda r: r.id, selected_seq_records))

        seqIDs_with_long_repeats = set(
            filter(
                lambda r: r.id in seqIDs_with_large_repeats,
                selected_seq_records
            )
        )

        num_genes_in_genome = gene_count_dict[ass_id]
        category = cat_dict[ass_id]

        if (num_genes_in_genome - len(seqIDs_with_long_repeats) == 0) and (category != 3):
            seqIDs_to_output = seqIDs
        else:
            seqIDs_to_output = seqIDs - seqIDs_with_long_repeats
        # end if

        seq_records_to_output = tuple(
            map(
                lambda r: r.id in seqIDs_to_output:
                selected_seq_records
            )
        )
        for record in seq_records_to_output:
            fasta_outfile.write(f'>{record.description}\n{str(record.seq)}\n')
        # end for
    # end for
# end with

os.system(f'seqkit stats -a {pure_genes_fpath}')

print('Calculating statistics')
gene_seqs_2_stats(pure_genes_fpath, ass_acc_fpath, pure_genes_stats_fpath)
print()

print('\nCompleted!')
