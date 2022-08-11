#!/usr/bin/env python3

# The script takes output file of the script `get_taxIDs.py` and output of the script
#   `pergenome_2_pergene_taxIDs.py` and uses them to map Assembly IDs and seqIDs to taxonomy
#   from file `rankedlineage.dmp` (argument `--ranked-lineage`) from,
#   see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump .
# Taxonomy is taxon names at the following levels:
#   superkingdom (domain), phylum, class, order, family, genus, species;
#   also term "taxonomy" here include taxIDs and NCBI taxonomy names.

## Command line arguments

### Input files:
# 1. `-i / --assm-acc-file` -- a TSV file (with header)
#   with Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs.
#   Mandatory.
# 1. `--per-genome-taxid-file` -- a TSV file mapping Assembly IDs to taxIDs.
#   This file is the output of the script `get_taxIDs.py`. Mandatory.
# 2. `--per-gene-taxid-file` -- a TSV file mapping genes seqIDs to taxIDs.
#   This file is the output of the script `pergenome_2_pergene_taxIDs.py`. Mandatory.
# 3. `--ranked-lineage` -- the file `rankedlineage.tsv` from NCBI's new_taxdump archive. Mandatory.

### Output files:
# 1. `--per-genome-outfile` -- a TSV file mapping Assembly IDs to taxonomy. Mandatory.
# 2. `--per-gene-outfile` -- a TSV file mapping RiboGrove seqIDs to taxonomy. Mandatory.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import re
import sys
import time
import socket
import argparse
import http.client

import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

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

args = parser.parse_args()


# For convenience
ass_acc_fpath = os.path.abspath(args.assm_acc_file)
per_genome_taxid_fpath = os.path.abspath(args.per_genome_taxid_file)
per_gene_taxid_fpath = os.path.abspath(args.per_gene_taxid_file)
rankedlineage_path = os.path.abspath(args.ranked_lineage)
per_genome_outfpath = os.path.abspath(args.per_genome_outfile)
per_gene_outfpath = os.path.abspath(args.per_gene_outfile)


# Check existance of all input files and dependencies
for fpath in (per_genome_taxid_fpath, per_gene_taxid_fpath, rankedlineage_path):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

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


def lingering_https_get_request(server, url, request_for=None, acc=None):
    # Function performs a "lingering" HTTPS request.
    # It means that the function tries to get the response
    #     again and again if the request fails.
    #
    # :param server: server address;
    # :type server: str;
    # :param url: the rest of url;
    # :type url: str;
    # :param request_for: some comment for error message;
    # :type request_for: str;
    # :param acc: target accession or some ID;
    # :type acc: str;
    #
    # Returns obtained response coded in UTF-8 ('str').

    error = True

    # We can get spurious 404 or sth due to instability of NCBI servers work.
    # Let's give it 3 attempts (with 15 sec spans in between),
    #   and if all them are unsuccessful -- teminate execution.
    attempt_i = 0
    max_attempts = 3

    while error:
        try:
            conn = http.client.HTTPSConnection(server, timeout=30) # create connection
            conn.request('GET', url) # ask for if there areresults
            response = conn.getresponse() # get the resonse

            if response.code != 200:
                if attempt_i < max_attempts and 'ncbi.nlm.nih.gov' in server:
                    print('Error {}: {}.'.format(response.code, response.reason))
                    print('It may be due to instable work of NCBI servers.')
                    print('{} attempts to connect left, waiting 15 sec...'\
                        .format(max_attempts - attempt_i))
                    attempt_i += 1
                else:
                    print('Cannot find {} for {}.'.format(request_for, acc))
                    print('Request failed with status code {}: {}'\
                        .format(response.code, response.reason))
                    return ''
                # end if
            # end if

            resp_content = str(response.read(), 'utf-8') # get response text
        except (OSError,\
                http.client.RemoteDisconnected,\
                socket.gaierror,\
                http.client.CannotSendRequest) as err:
            comment_str = ''
            if not request_for is None:
                comment_str += ' requesting for {}'.format(request_for)
                if not acc is None:
                    comment_str += ' (accession: `{}`)'.format(acc)
                # end if
                comment_str += '.'
            # end if
            print()
            print('Can\'t connect to `{}`{}'.format(server + url, comment_str))
            print( str(err) )
            print('the program will sleep for 30 seconds and try to connect again.')
            time.sleep(30)
        else:
            error = False # if no exception ocured, get out of the loop
        finally:
            conn.close()
        # end try
    # end while
    return resp_content
# end def lingering_https_get_request


def request_missing_taxonomy(taxid):
    # Function retrieves taxonomy of a hit from NCBI.
    #
    # :param taxid: target Taxonomy ID;
    # :type taxid: int;

    ranks = ('superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')

    # These words at second (with index 1) position of title indicate that
    #   actual species name are specified after it.
    second_words_not_species = ('species', 'sp.', 'strain', 'str.', 'bacterium')

    # Get taxonomy page of the organism
    taxonomy_url = '/Taxonomy/Browser/wwwtax.cgi?mode=Info&id={}'.format(taxid)
    taxonomy_text = lingering_https_get_request('www.ncbi.nlm.nih.gov',
        taxonomy_url, 'taxonomy', taxid)

    # This pattern will match taxonomic names along with their ranks
    tax_rank_pattern = r'TITLE=\"([a-z ]+)\"\>([A-Z].+?)\</a\>'

    # Get all taxonomic names of the organism
    taxonomy = re.findall(tax_rank_pattern, taxonomy_text)

    # E.g., this record has no appropriate ranks: CP034535
    # Merely return it's definition
    if len(taxonomy) == 0:
        return {
            'tax_name': 'NA',
            'domain': 'NA',
            'phylum': 'NA',
            'class': 'NA',
            'order': 'NA',
            'family': 'NA',
            'genus': 'NA',
            'species': 'NA',
        }
    # end if

    # We will convert ranks to lowercase just in case.
    taxonomy_dict = {
        tax_tuple[0].lower(): tax_tuple[1] for tax_tuple in taxonomy
    }

    # We will leave only following taxonomic ranks: domain, phylum, class, order, family, genus.
    # Species name requires special handling, it will be added later.
    ranks_to_select = ranks[:-1]

    # Remove redundant ranks
    all_available_ranks = set(taxonomy_dict.keys())
    for rank_name in all_available_ranks:
        if not rank_name in ranks_to_select:
            del taxonomy_dict[rank_name]
        # end if
    # end for

    # Parse and add taxonomy name to the dictionary
    tax_name = re.search(r"\<title\>Taxonomy browser \((.+)\)\</title\>", taxonomy_text).group(1)
    taxonomy_dict['tax_name'] = tax_name

    # Check if species name is specified like other ranks:
    check_direct_species_patt = r'TITLE=\"(species)\"\>([A-Za-z0-9 \.]+)\</a\>'
    match_direct_species = re.search(check_direct_species_patt, taxonomy_text)

    if not match_direct_species is None:
        # If species name is specified like other ranks, merely add it to list:
        taxonomy_dict['species'] = taxonomy_dict['genus'] + ' ' + match_direct_species.group(2).partition(' ')[2]
    else:
        # Otherwise we need to parse species name from tax_name

        # Get words
        tax_name_words = tax_name.split(' ')

        # We will take all this words as species name.
        #   Example: MN908947
        try:
            if tax_name_words[1] in second_words_not_species:
                taxonomy_dict['species'] = taxonomy_dict['genus'] + ' ' + ' '.join(tax_name_words[1:])
            else:
                taxonomy_dict['species'] = taxonomy_dict['genus'] + ' ' + tax_name_words[1]
            # end if
        except IndexError:
            # Handle absence of species name, e.g., this: AC150248.3
            # Well, nothing to append in this case!
            pass
        # end try
    # end if

    # Fill in missing ranks with NAs
    obtained_ranks = set(taxonomy_dict.keys())
    for rank_name in ranks:
        if not rank_name in obtained_ranks:
            taxonomy_dict[rank_name] = 'NA'
        # end if
    # end for

    # Rename 'superkingdom' -> 'domain'
    taxonomy_dict['domain'] = taxonomy_dict['superkingdom']
    del taxonomy_dict['superkingdom']

    return taxonomy_dict
# end def request_missing_taxonomy


def fill_missing_taxonomy(row):

    if pd.isnull(row['tax_name']):

        print(f'Requesting taxonomy for taxID {row["taxID"]}... ')

        taxonomy_dict = request_missing_taxonomy(row['taxID'])

        # If no taxonomy was retrieved, get it from the RefSeq title
        # Bad way, but no better ways are left
        if taxonomy_dict['tax_name'] == 'NA':
            global ass_acc_df
            seq_title = ass_acc_df[ass_acc_df['ass_id'] == row['ass_id']] \
                .reset_index().loc[0, 'title']
            strings_to_rm = (
                ', complete sequence',
                ', complete genome',
                ' map unlocalized',
            )
            for str_to_rm in strings_to_rm:
                seq_title = seq_title.replace(str_to_rm, '')
            # end for
            taxonomy_dict['tax_name'] = seq_title

            print(f'Cannot find taxonomy for taxID {row["taxID"]} at the NCBI website')
            print(f'Using the RefSeq title as the taxonomy name: `{taxonomy_dict["tax_name"]}`')
        else:
            print(taxonomy_dict)
        # end if

        # Fill the taxonomy
        for rank_name, taxon_name in taxonomy_dict.items():
            row[rank_name] = taxon_name
        # end for
    # end if

    return row
# end def fill_missing_taxonomy


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

# Rename the 'superkingdom' column
rankedlineage_df = rankedlineage_df.rename(columns={'superkingdom': 'domain'})

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


# Sometimes information is missing for some taxIDs in rankedlineade.dmp
# Request the missing taxonomy from NCBI Taxonomy
missing_taxIDs = set(
    per_genome_taxonomy_df[
        pd.isnull(per_genome_taxonomy_df['tax_name'])
    ]['taxID']
)
if len(missing_taxIDs) != 0:
    print(f'Taxonomy is missing for {len(missing_taxIDs)} Taxonomy IDs')
    print('The script will request the taxonomy for them from the NCBI website')
    ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t')
    per_genome_taxonomy_df = per_genome_taxonomy_df.apply(fill_missing_taxonomy, axis=1)
    del ass_acc_df
# end if
del missing_taxIDs


# Amend class for Cyanobacteria
per_genome_taxonomy_df = per_genome_taxonomy_df.apply(amend_Cyanophyceae, axis=1)
# Amend species names
per_genome_taxonomy_df = per_genome_taxonomy_df.apply(fill_empty_species_name, axis=1)

# Order columns
per_genome_taxonomy_df = per_genome_taxonomy_df[
    [
        'ass_id',
        'accs',
        'taxID',
        'tax_name',
        'species',
        'genus',
        'family',
        'order',
        'class',
        'phylum',
        'domain',
    ]
]

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
per_genome_df_for_merging = per_genome_taxonomy_df[
    [
        'taxID', 'tax_name',
        'domain', 'phylum', 'class', 'order', 'family', 'genus', 'species',
    ]
].drop_duplicates()

per_gene_taxonomy_df = per_gene_taxid_df.merge(
    # rankedlineage_df,
    per_genome_df_for_merging,
    on='taxID',
    how='left'
)
# Amend class for Cyanobacteria
per_gene_taxonomy_df = per_gene_taxonomy_df.apply(amend_Cyanophyceae, axis=1)
# Amend species names
per_gene_taxonomy_df = per_gene_taxonomy_df.apply(fill_empty_species_name, axis=1)

# Order columns
per_gene_taxonomy_df = per_gene_taxonomy_df[
    [
        'ass_id',
        'seqID',
        'taxID',
        'tax_name',
        'species',
        'genus',
        'family',
        'order',
        'class',
        'phylum',
        'domain',
    ]
]


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
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
