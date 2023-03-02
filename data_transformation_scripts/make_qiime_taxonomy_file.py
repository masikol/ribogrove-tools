#!/usr/bin/env python3

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import re
import sys
import argparse

import numpy as np
import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--in-ribogrove-gene-stats',
    help='RiboGrove per-gene statistics tsv file `gene_seqs_statistics.tsv`',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--out-qiime2-taxonomy',
    help='output file -- QIIME2-compatible taxonomy',
    required=True
)

args = parser.parse_args()


# For convenience
in_gene_stats_fpath = os.path.abspath(args.in_ribogrove_gene_stats)
outfpath = os.path.abspath(args.out_qiime2_taxonomy)


# Check existance of all input files and dependencies
if not os.path.exists(in_gene_stats_fpath):
    print(f'Error: file `{in_gene_stats_fpath}` does not exist!')
    sys.exit(1)
# end if

# Create output directory if needed
if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

print(in_gene_stats_fpath)
print()

def make_qiime2_outrow(rg_row):
    out_dict = {
        'seqID': rg_row['seqID'],
        'taxonomy': make_qiime2_tax_string(rg_row),
    }
    return pd.Series(out_dict)
# end def

def make_qiime2_tax_string(rg_row):
    # k__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; o__Legionellales; f__Legionellaceae; g__Legionella; s__
    qiime2_tax_str = 'k__{}; p__{}; c__{}; o__{}; f__{}; g__{}; s__{}'.format(
        '' if pd.isnull(rg_row['Domain']) else rg_row['Domain'],
        '' if pd.isnull(rg_row['Phylum']) else rg_row['Phylum'],
        '' if pd.isnull(rg_row['Class'])  else rg_row['Class'],
        '' if pd.isnull(rg_row['Order'])  else rg_row['Order'],
        '' if pd.isnull(rg_row['Family']) else rg_row['Family'],
        '' if pd.isnull(rg_row['Genus'])  else rg_row['Genus'],
        make_species_name(
            '' if pd.isnull(rg_row['Species']) else rg_row['Species']
        )
    )
    return qiime2_tax_str
# end def

def make_species_name(raw_species):
    species_name = raw_species
    
    if raw_species.startswith('Candidatus '):
        species_name = species_name.replace('Candidatus ', '')
    # end if

    species_name = species_name.partition(' ')[2]

    return species_name
# end def


rg_gene_stats_df = pd.read_csv(in_gene_stats_fpath, sep='\t')

out_df = pd.DataFrame(
    {
        'seqID': np.repeat('', rg_gene_stats_df.shape[0]),
        'taxonomy': np.repeat('', rg_gene_stats_df.shape[0]),
    }
)

i = 0
inc_i = rg_gene_stats_df.shape[0] // 100
next_i = inc_i

sys.stdout.write('0/{} done, 0%'.format(rg_gene_stats_df.shape[0]))
sys.stdout.flush()

for i, row in rg_gene_stats_df.iterrows():
    out_series = make_qiime2_outrow(row)
    out_df.iloc[i,] = out_series

    i += 1
    if i >= next_i:
        sys.stdout.write(
            '\r{}/{} done, {}%'.format(
                i,
                rg_gene_stats_df.shape[0],
                round(i / rg_gene_stats_df.shape[0] * 100, 2)
            )
        )
        sys.stdout.flush()
        next_i += inc_i
    # end if
# end for

print(
    '\r{}/{} done, {}%'.format(
        i,
        rg_gene_stats_df.shape[0],
        round(i / rg_gene_stats_df.shape[0] * 100, 2)
    )
)

out_df.to_csv(
    outfpath,
    sep='\t',
    na_rep='',
    header=False,
    index=False,
    encoding='utf-8'
)


print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
