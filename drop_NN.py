#!/usr/bin/env python3

import re

from Bio import SeqIO

from gene_seqs_2_stats import gene_seqs_2_stats


# seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected.fasta'
seqs_fpath = '/home/deynonih/cager/new_16S_scrubbling/gene_seqs/all_collected.fasta'

# ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
ass_acc_fpath = '/home/deynonih/cager/new_16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'

# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'
outfpath = '/home/deynonih/cager/new_16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'

# nn_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/NN/gene_seqs_NN.fasta'
nn_seqs_fpath = '/home/deynonih/cager/new_16S_scrubbling/gene_seqs/gene_seqs_NN.fasta'

# out_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/out_stats_no_NN.tsv'
out_stats_fpath = '/home/deynonih/cager/new_16S_scrubbling/gene_seqs/out_stats_no_NN.tsv'


nn_pattern = r'NN'
nn_count = 0

next_report = 499
inc = 500


num_seqs = len(tuple(SeqIO.parse(seqs_fpath, 'fasta')))
seq_records = SeqIO.parse(seqs_fpath, 'fasta')


with open(outfpath, 'wt') as outfile, open(nn_seqs_fpath, 'wt') as outfile_nn:

    for i, record in enumerate(seq_records):

        if i == next_report:
            print(f'\r{i+1}/{num_seqs}', end=' ')
            next_report += inc
        # end if

        if re.search(nn_pattern, str(record.seq)) is None:
            outfile.write(f'>{record.id}\n{str(record.seq)}\n')
        else:
            nn_count += 1
            outfile_nn.write(f'>{record.id}\n{str(record.seq)}\n')
        # end if
    # end for
# end with

print(f'\r{i+1}/{num_seqs}\n')
print(f'{nn_count} NN-sequences found')
print(outfpath)
print(nn_seqs_fpath)

print('Calculating statistics')
gene_seqs_2_stats(outfpath, ass_acc_fpath, out_stats_fpath)
print(out_stats_fpath)

print('Completed!')
