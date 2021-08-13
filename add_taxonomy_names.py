#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import numpy as np
import pandas as pd


infpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/taxIDs.tsv'
per_gene_infpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/per_gene_taxIDs.tsv'

rankedlineage_path = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/new-taxdump/rankedlineage_2.dmp'

outfpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/taxonomy.tsv'
per_gene_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/per_gene_taxonomy.tsv'



def amend_Cyanophyceae(row):

    if row['phylum'] == 'Cyanobacteria':
        row['class'] = 'Cyanophyceae'
    # end if

    return row
# end def amend_Cyanophyceae


# =========

rankedlineage_df = pd.read_csv(
    rankedlineage_path,
    sep='\t',
    names=['taxID', 'tax_name', 'species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom', 'superkingdom'],
    header=None,
    dtype={
        'taxID': pd.Int32Dtype(),
        'tax_name': str,
        'species': str,
        'genus': str,
        'family': str,
        'order': str,
        'class': str,
        'phylum': str,
        'kingdom': str,
        'superkingdom': str
    }
)

rankedlineage_df = rankedlineage_df.drop(columns=['species', 'kingdom'], axis=1)

# =========

print('Merging input file without seqIDs')

taxid_df = pd.read_csv(
    infpath,
    sep='\t',
    dtype={
        'ass_id': int,
        'accs': str,
        'taxID': pd.Int32Dtype()
    }
)


taxonomy_df = taxid_df.merge(rankedlineage_df, on='taxID', how='left')

taxonomy_df = taxonomy_df.apply(amend_Cyanophyceae, axis=1)

print(taxonomy_df.shape)
print(taxonomy_df.head())

taxonomy_df.to_csv(
    outfpath,
    sep='\t',
    na_rep='NA',
    header=True,
    index=False,
    encoding='utf-8'
)

print(outfpath)

# =========

print('Merging input file with seqIDs..')

per_gene_taxid_df = pd.read_csv(
    per_gene_infpath,
    sep='\t',
    dtype={
        'seqID': str,
        'ass_id': int,
        'accs': str,
        'taxID': pd.Int32Dtype()
    }
)

per_gene_taxonomy_df = per_gene_taxid_df.merge(rankedlineage_df, on='taxID', how='left')

per_gene_taxonomy_df = per_gene_taxonomy_df.apply(amend_Cyanophyceae, axis=1)

print(per_gene_taxonomy_df.shape)
print(per_gene_taxonomy_df.head())

per_gene_taxonomy_df.to_csv(
    per_gene_outfpath,
    sep='\t',
    na_rep='NA',
    header=True,
    index=False,
    encoding='utf-8'
)
print(per_gene_outfpath)

print('Completed!')
