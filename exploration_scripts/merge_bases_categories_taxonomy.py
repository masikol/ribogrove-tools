#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import pandas as pd

# Bacteria
bases_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bases_count.tsv'
categories_fpath = '/mnt/1.5_drive_0/16S_scrubbling/categories/bacteria_per_gene_categories.tsv'
taxonomy_fpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/per_gene_taxonomy.tsv'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/pure_genes_per_gene_stats.tsv'

# Archaea
# bases_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bases_count.tsv'
# categories_fpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/categories/archaea_per_gene_categories.tsv'
# taxonomy_fpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/per_gene_taxIDs.tsv'
# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/pure_genes_per_gene_stats.tsv'


# seqID   a       t       g       c       len
bases_df = pd.read_csv(bases_fpath, sep='\t')

# ass_id  seqID   category
categories_df = pd.read_csv(categories_fpath, sep='\t')

# seqID   ass_id  accs    taxID   tax_name        genus   family  order   class   phylum  superkingdom
taxonomy_df = pd.read_csv(taxonomy_fpath, sep='\t')


print('Merging...')
merged_df = bases_df \
    .merge(categories_df, on='seqID', how='left') \
    .merge(taxonomy_df.drop(['ass_id'], axis=1), on='seqID', how='left')
print('  done!\n')

print('Merged dataframe:')
print(merged_df.shape)
print(merged_df.head())

merged_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    encoding='utf-8',
    na_rep='NA'
)

print('Completed!')
print(outfpath)
