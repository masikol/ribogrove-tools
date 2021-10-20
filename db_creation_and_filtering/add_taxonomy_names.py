#!/usr/bin/env python3

# The script takes both output files of the script `get_taxIDs.py` (--per-genome-taxid-file, --per-gene-taxid-file)
#   and uses them to map Assemblty IDs and seqIDs to taxonomy from file rankedlineage.dmp (--ranked-lineage).

# Input files:
# 1. TSV file mapping Assembly IDs to taxIDs (--per-genome-taxid-file)
# 2. TSV file mapping genes seqIDs to taxIDs (--per-gene-taxid-file)
# 3. File `rankedlineage.tsv` from NCBI's new_taxdump archive
#   Option --ranked-lineage.
#   See https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/

# Output files:
# 1. TSV file mapping Assembly IDs to taxonomy (--per-genome-outfile)
# 2. TSV file mapping genes seqIDs to taxonomy (--per-gene-outfile)

# Dependencies:
# 1. Path to seqkit exetutable (--seqkit)


import os
import sys
import argparse

import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '--per-genome-taxid-file',
    help="""TSV file (with header) mapping Assembly IDs to taxIDs.
    It should have 3 columns (ass_id, accs, taxID).""",
    required=True
)

parser.add_argument(
    '--per-gene-taxid-file',
    help="""TSV file (with header) mapping genes SeqIDs to taxIDs.
    It should have 4 columns (seqID, ass_id, accs, taxID).""",
    required=True
)

parser.add_argument(
    '--ranked-lineage',
    help="""file rankedlineage.tsv from NCBI\'s new_taxdump archive
    see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz""",
    required=True
)

# Output files

parser.add_argument(
    '--per-genome-outfile',
    help='output file mapping Assembly IDs to taxonomy',
    required=True
)

parser.add_argument(
    '--per-gene-outfile',
    help='output file mapping genes seqIDs to taxonomy',
    required=True
)

# Dependencies

parser.add_argument(
    '--seqkit',
    help='seqkit executable',
    required=True
)

args = parser.parse_args()


# For convenience
per_genome_taxid_fpath = os.path.abspath(args.per_genome_taxid_file)
per_gene_taxid_fpath = os.path.abspath(args.per_gene_taxid_file)
rankedlineage_path = os.path.abspath(args.ranked_lineage)
per_genome_outfpath = os.path.abspath(args.per_genome_outfile)
per_gene_outfpath = os.path.abspath(args.per_gene_outfile)
seqkit_fpath = os.path.abspath(args.seqkit)


# Check existance of all input files and dependencies
for fpath in (per_genome_taxid_fpath, per_gene_taxid_fpath, rankedlineage_path, seqkit_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Check if seqkit executable is actually executable
if not os.access(seqkit_fpath, os.X_OK):
    print(f'Error: file `{seqkit_fpath}` is not executable!')
    sys.exit(1)
# end if

# Create output directories if needed
for some_dir in map(os.path.dirname, [per_genome_outfpath, per_gene_outfpath]):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end if

print(per_genome_taxid_fpath)
print(per_gene_taxid_fpath)
print(rankedlineage_path)
print(seqkit_fpath)
print()


def reformat_rankedlineaage_file(rankedlineage_path: str) -> str:
    # Original rankedlineage.dmp file has very weird separator: "\t|\t". And terminal "\t|" (OMG why???).
    # We will replace "\t|\t" with mere "\t" and remove terminal "\t|".

    new_rankedlineade_fpath = os.path.join(
        os.path.dirname(rankedlineage_path),
        'rankedlineage_just_tabs.dmp'
    )

    cmd = f'cat {rankedlineage_path} | sed "s/\\t|\\t/\\t/g" | sed "s/\\t|//g" > {new_rankedlineade_fpath}'
    print(cmd)
    exit_code = os.system(cmd)

    if exit_code != 0:
        print(f'Error: cannot reformat rankedlineage file: {rankedlineage_path}')
        sys.exit(1)
    # end if

    return new_rankedlineade_fpath
# end def reformat_rankedlineaage_file



def amend_Cyanophyceae(row: pd.core.series.Series) -> pd.core.series.Series:
    # NCBI taxonomy violently omits 'class' taxon for Cyanobacteria.
    # But they have class: all the same -- Cyanophyceae (see https://lpsn.dsmz.de/class/cyanophyceae).
    # This function restores justice. Long live Cyanobacteria!
    # The function is meant to be used with `pd.DataFrame.apply` function

    if row['phylum'] == 'Cyanobacteria':
        row['class'] = 'Cyanophyceae'
    # end if

    return row
# end def amend_Cyanophyceae


def fill_empty_species_name(row: pd.core.series.Series) -> pd.core.series.Series:
    # If taxID points to a species, (like 1642), it's taxonomy contains no 'species' fiels.
    # Well, then we will copy it from 'tax_name' field.

    if pd.isnull(row['species']):
        row['species'] = row['tax_name']
    # end if

    return row
# end def fill_empty_species_name


# == Proceed ==

# Reformat rankedlineage file
print(f'Reformatting file `{rankedlineage_path}` in order to make it manageable...')
reformatted_rankedlineade_fpath = reformat_rankedlineaage_file(rankedlineage_path)
print('done.\n')


# Read reformatted rankedlineage file

print('Reading reformatted rankedlineage file...')

rankedlineage_df = pd.read_csv(
    reformatted_rankedlineade_fpath,
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

# Remove columns of no interest
rankedlineage_df = rankedlineage_df.drop(columns=['kingdom'], axis=1)


# Make per-genome taxonomy file

print('Creating per-genome taxonomy file')

# Read per-genome taxID file
taxid_df = pd.read_csv(
    per_genome_taxid_fpath,
    sep='\t',
    dtype={
        'ass_id': int,
        'accs': str,
        'taxID': pd.Int32Dtype()
    }
)

# Merge per-genome taxID file to rankedlineage file
per_genome_taxonomy_df = taxid_df.merge(rankedlineage_df, on='taxID', how='left')
# Amend class for Cyanobacteria
per_genome_taxonomy_df = per_genome_taxonomy_df.apply(amend_Cyanophyceae, axis=1)
# Amend species names
per_genome_taxonomy_df = per_genome_taxonomy_df.apply(fill_empty_species_name, axis=1)

# Write output per-genome file
per_genome_taxonomy_df.to_csv(
    per_genome_outfpath,
    sep='\t',
    na_rep='NA',
    header=True,
    index=False,
    encoding='utf-8'
)


# Make per-gene taxonomy file

print('Creating per-gene taxonomy file')

# Read per-gene taxID file
per_gene_taxid_df = pd.read_csv(
    per_gene_taxid_fpath,
    sep='\t',
    dtype={
        'seqID': str,
        'ass_id': int,
        'accs': str,
        'taxID': pd.Int32Dtype()
    }
)

# Merge per-gene taxID file to rankedlineage file
per_gene_taxonomy_df = per_gene_taxid_df.merge(rankedlineage_df, on='taxID', how='left')
# Amend class for Cyanobacteria
per_gene_taxonomy_df = per_gene_taxonomy_df.apply(amend_Cyanophyceae, axis=1)
# Amend species names
per_gene_taxonomy_df = per_gene_taxonomy_df.apply(fill_empty_species_name, axis=1)

# Write output per-gene file
per_gene_taxonomy_df.to_csv(
    per_gene_outfpath,
    sep='\t',
    na_rep='NA',
    header=True,
    index=False,
    encoding='utf-8'
)

print('Completed!')
print(per_genome_outfpath)
print(per_gene_outfpath)
