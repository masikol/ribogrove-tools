#!/usr/bin/env python3

import re

# https://github.com/deprekate/RepeatFinder
import repeatfinder as rf
from Bio import SeqIO


seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/repeats_no_NN.tsv'


nonredundant_pattern = r'[ATGC]'
next_report = 499
inc = 500

num_seqs = len(tuple(SeqIO.parse(seqs_fpath, 'fasta')))
seq_records = SeqIO.parse(seqs_fpath, 'fasta')


def get_repeat_len(repeat_out):
    return repeat_out[1] - repeat_out[0] + 1
# end def get_repeat_len

with open(outfpath, 'wt') as outfile:

    outfile.write(f'seqID\tr1_start\tr1_end\tr2_start\tr2_end\trep_len\n')

    for i, record in enumerate(seq_records):

        if i == next_report:
            print(f'\r{i+1}/{num_seqs}', end=' ')
            next_report += inc
        # end if

        repeats = rf.get_repeats(str(record.seq))

        for r in repeats:
            if not re.search(nonredundant_pattern, str(record.seq)[r[0]-1 : r[1]]) is None:
                rep_len = get_repeat_len(r)
                outfile.write(f'{record.id}\t{r[0]}\t{r[1]}\t{r[2]}\t{r[3]}\t{rep_len}\n')
            # end if
        # end for
    # end for
# end with

print(f'\r{i+1}/{num_seqs}')
print('Completed!')
print(outfpath)
