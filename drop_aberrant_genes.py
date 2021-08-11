#!/usr/bin/env python3

import os


all_genes_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected.fasta'
aberrant_seqIDs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/aberrant_seqIDs.txt'

no_aberrants_genes_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/pure_genes.fasta'

os.system(f'seqkit stats {all_genes_seqs_fpath}')

cmd = f'cat {all_genes_seqs_fpath} | seqkit grep -vf {no_aberrants_genes_seqs_fpath} > {no_aberrants_genes_seqs_fpath}'
print(cmd)
os.system(cmd)

os.system(f'seqkit stats {no_aberrants_genes_seqs_fpath}')

print('\nCompleted!')
print(no_aberrants_genes_seqs_fpath)
