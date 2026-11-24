#!/usr/bin/env python3

# The script makes a taxonomy file for the downloaded genomes.

## Command line arguments

### Input files:
# 1. `-i / --asm-sum` -- an assembly summary file after the 2nd step of filtering.
#   Mandatory.
# 2. `-l / --ranked-lineage` -- file rankedlineage.tsv from NCBI's new_taxdump archive
#   see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
#   Mandatory.
# 3. `-d / --domain` -- domain name: "bacteria" or "archaea".
#   Mandatory.


### Output files:
# 1. `-o / --out` -- output taxonomy file.
#   Mandatory.


import os
from src.rg_tools_time import get_time

print(
    '\n|=== {} STARTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)


# == Parse arguments ==
import argparse

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--asm-sum',
    help='an assembly summary file after the 2nd step of filtering',
    required=True
)

parser.add_argument(
    '-l',
    '--ranked-lineage',
    help="""file rankedlineage.tsv from NCBI\'s new_taxdump archive
    see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz""",
    required=True
)

parser.add_argument(
    '-d',
    '--domain',
    help='domain name, for example, "bacteria"',
    required=True
)

# Output file

parser.add_argument(
    '-o',
    '--out',
    help='output taxonomy file',
    required=True
)


args = parser.parse_args()


# == Import them now ==
import re
import sys
import http
import time
import socket

import polars as pl
from Bio import Entrez

import src.rg_tools_IO as rgIO

# NCBI requires an email for Entrez requests
Entrez.email = 'sikolenko@mbio.bas-net.by'

# For convenience
asm_sum_fpath = os.path.abspath(args.asm_sum)
rankedlineage_path = os.path.abspath(args.ranked_lineage)
domain = args.domain.capitalize()
outfpath = os.path.abspath(args.out)


# Check existance of all input files and dependencies
for fpath in (asm_sum_fpath, rankedlineage_path):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for


# Create output directories if needed
if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        print(err)
        sys.exit(1)
    # end try
# end if


def reformat_rankedlineaage_file(rankedlineage_path: str) -> str:
    # Original rankedlineage.dmp file has very weird separator: "\t|\t". And terminal "\t|" (OMG why???).
    # We will replace "\t|\t" with mere "\t" and remove terminal "\t|".

    new_rankedlineage_fpath = os.path.join(
        os.path.dirname(rankedlineage_path),
        'rankedlineage_just_tabs.dmp'
    )

    cmd = ' '.join([
        'cat', rankedlineage_path,
        '|',
        'sed "s/\\t|\\t/\\t/g"',
        '|',
        'sed "s/\\t|//g"',
        '|',
        'tr -d \'"\'',
        '>', new_rankedlineage_fpath
    ])
    print(cmd)
    exit_code = os.system(cmd)

    if exit_code != 0:
        print(f'Error: cannot reformat rankedlineage file: {rankedlineage_path}')
        sys.exit(1)
    # end if

    return new_rankedlineage_fpath
# end def

def fill_empty_species_name(row):
    # If taxid points to a species, (like 1642), it's taxonomy contains no 'species' fiels.
    # Well, then we will copy it from 'organism_name' field.

    if row['Species'] is None:
        row['Species'] = row['organism_name']
    # end if

    return row['Species']
# end def


def get_lineage_for_taxid(taxid: str) -> dict:
    ranks_of_interest = (
        'Domain',
        'Kingdom',
        'Phylum',
        'Class',
        'Order',
        'Family',
        'Genus',
        'Species'
    )

    result = {rank: None for rank in ranks_of_interest}

    # Fetch the taxonomy record as XML
    error = True

    # We can get spurious 404 or sth due to instability of NCBI servers work.
    # Let's give it 3 attempts (with 15 sec spans in between),
    #   and if all them are unsuccessful -- teminate execution.
    attempt_i = 0
    max_attempts = 3

    while error:
        try:
            handle = Entrez.efetch(db='taxonomy', id=taxid, retmode='xml')
            records = Entrez.read(handle)
        except (OSError,\
                http.client.RemoteDisconnected,\
                socket.gaierror,\
                http.client.CannotSendRequest) as err:
            print()
            print("Can't connect to Entrez requestring for taxid {}".format(taxid))
            print(str(err))
            print('The program will sleep for 30 seconds and try to connect again.')
            time.sleep(30)
        else:
            error = False # if no exception ocured, get out of the loop
        finally:
            handle.close()
        # end try
    # 3nd while

    if not records:
        return result
    # end if

    record = records[0]
    result['organism_name'] = record['ScientificName']

    # NCBI uses 'superkingdom' for what is strictly called Domain
    rank_aliases = {
        'superkingdom': 'Domain',
        'domain': 'Domain',
    }

    # LineageEx is a list of dicts, each with 'Rank' and 'ScientificName'
    for ancestor in record.get('LineageEx', []):
        ncbi_rank = ancestor.get('Rank', '').lower()
        name = ancestor.get('ScientificName', '')

        # Map NCBI rank to our requested rank names
        mapped_rank = rank_aliases.get(ncbi_rank, ncbi_rank.capitalize())

        if mapped_rank in result and result[mapped_rank] is None:
            result[mapped_rank] = name
        # end if
    # end for

    # The query organism itself is not in LineageEx; add it if it matches a requested rank
    query_rank = record.get('Rank', '').lower()
    mapped_query_rank = rank_aliases.get(query_rank, query_rank.capitalize())
    if mapped_query_rank in result and result[mapped_query_rank] is None:
        result[mapped_query_rank] = record.get('ScientificName')
    # end if

    return result
# end def


def fill_missing_taxonomy(row: tuple) -> tuple:
    # Column order:
    # asm_acc
    # taxid
    # organism_name
    # Species
    # Genus
    # Family
    # Order
    # Class
    # Phylum
    # Kingdom
    # Domain
    organism_name = row[2]

    if organism_name is None:
        asm_acc = row[0]
        taxid = row[1]

        print(f'Requesting taxonomy for taxid {taxid}... ')
        taxonomy_dict = get_lineage_for_taxid(taxid)
        
        if taxonomy_dict['organism_name'] == 'NA':
            raw_org_name = asm_sum_df.filter(pl.col('asm_acc') == asm_acc)['raw_organism_name'][0]
            strings_to_rm = (
                ', complete sequence',
                ', complete genome',
                ' map unlocalized',
            )
            for str_to_rm in strings_to_rm:
                raw_org_name = raw_org_name.replace(str_to_rm, '')
            # end for

            taxonomy_dict['organism_name'] = raw_org_name
            print(f'Cannot find taxonomy for taxid {taxid} at the NCBI website')
            print(f'Using the RefSeq title as the taxonomy name: `{taxonomy_dict["organism_name"]}`')
            print('The script will terminate with a non-ero status code.')
            print('Please fill the missing taxonomy data manually')
            global EXIT_CODE
            EXIT_CODE = 1
        else:
            print(taxonomy_dict)
        # end if

        # Update row with taxonomy values
        for rank_name, taxon_name in taxonomy_dict.items():
            taxonomy_dict[rank_name] = taxon_name
        # end for

        taxonomy_tuple = (
            asm_acc,
            taxid,
            taxonomy_dict['organism_name'],
            taxonomy_dict['Species'],
            taxonomy_dict['Genus'],
            taxonomy_dict['Family'],
            taxonomy_dict['Order'],
            taxonomy_dict['Class'],
            taxonomy_dict['Phylum'],
            taxonomy_dict['Kingdom'],
            taxonomy_dict['Domain'],
        )
        return taxonomy_tuple
    # end if
    return row
# end def


# == Proceed ==

EXIT_CODE = 0

# Reformat rankedlineage file
print(f'Reformatting file `{rankedlineage_path}` in order to make it manageable...')
reformatted_rankedlineage_fpath = reformat_rankedlineaage_file(rankedlineage_path)
print('done.\n')


# Read reformatted rankedlineage file

print('Reading reformatted rankedlineage file...')

rankedlineage_df = pl.read_csv(
    reformatted_rankedlineage_fpath,
    separator='\t',
    n_threads=1,
    new_columns=[
        'taxid', 'organism_name',
        'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum',
        'Kingdom', 'Domain'
    ],
    has_header=False,
    schema={
        'taxid': pl.UInt32,
        'organism_name': pl.String,
        'Species': pl.String,
        'Genus': pl.String,
        'Family': pl.String,
        'Order': pl.String,
        'Class': pl.String,
        'Phylum': pl.String,
        'Kingdom': pl.String,
        'Domain': pl.String,
    },
    null_values=['NA', 'na', '']
)


# Make per-genome taxonomy file

print('Creating taxonomy file')

# Read per-genome taxid file
asm_sum_df = rgIO.read_ass_sum_file(asm_sum_fpath)
# Rename obsolete organism_name: it will be substituted by that from rankedlineage_df
asm_sum_df = asm_sum_df.rename({'organism_name': 'raw_organism_name'})

# Merge per-genome taxid file to rankedlineage file
taxonomy_df = asm_sum_df \
    .join(rankedlineage_df, on='taxid', how='left') \
    .select(
        pl.col(
            'asm_acc',
            'taxid',
            'organism_name',
            'Species',
            'Genus',
            'Family',
            'Order',
            'Class',
            'Phylum',
            'Kingdom',
            'Domain',
        )
    )
del rankedlineage_df

# Amend species names
taxonomy_df = taxonomy_df.with_columns(
    pl.struct(['Species', 'organism_name']).map_elements(
        fill_empty_species_name,
        return_dtype=pl.String
    ).alias('Species')
)

# Sometimes information is missing for some taxids in rankedlineage.dmp
# Request the missing taxonomy from NCBI Taxonomy
missing_taxids = frozenset(
    taxonomy_df.filter(
        pl.col('organism_name').is_null()
    )['taxid']
)
if len(missing_taxids) != 0:
    print(f'Taxonomy is missing for {len(missing_taxids)} Taxonomy IDs')
    print('The script will request the taxonomy for them from the NCBI website')
    column_names = list(taxonomy_df.columns)
    taxonomy_df = taxonomy_df.map_rows(fill_missing_taxonomy)
    taxonomy_df.columns = column_names
# end if
del missing_taxids

print(taxonomy_df.columns)

# Write output per-genome file
taxonomy_df.write_csv(
    outfpath,
    separator='\t',
    null_value='NA',
    include_header=True
)


print('\nCompleted!')
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)

sys.exit(EXIT_CODE)
