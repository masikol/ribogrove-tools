#!/usr/bin/env python3

# The script makes a taxonomy file for the downloaded genomes.

## Command line arguments

### Input files:
# 1. `-i / --asm-sum` -- an assembly summary file after the 2nd step of filtering.
#   Mandatory.
# 2. `-l / --ranked-lineage` -- file rankedlineage.tsv from NCBI's new_taxdump archive
#   see https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/new_taxdump/new_taxdump.tar.gz
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
import sys

import pandas as pd

import src.rg_tools_IO as rgIO


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
        sys.exit(1)
    # end try
# end if


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
# end def



def amend_Cyanophyceae(row: pd.Series) -> pd.Series:
    # NCBI taxonomy violently omits 'class' taxon for Cyanobacteria.
    # But they have class: all the same -- Cyanophyceae (see https://lpsn.dsmz.de/class/cyanophyceae).
    # This function restores justice. Long live Cyanobacteria!
    # The function is meant to be used with `pd.DataFrame.apply` function

    if row['Phylum'] == 'Cyanobacteria':
        row['Class'] = 'Cyanophyceae'
    # end if

    return row
# end def


def fill_empty_species_name(row: pd.Series) -> pd.Series:
    # If taxid points to a species, (like 1642), it's taxonomy contains no 'species' fiels.
    # Well, then we will copy it from 'organism_name' field.

    if pd.isnull(row['Species']):
        row['Species'] = row['organism_name']
    # end if

    return row
# end def


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
# end def


def request_missing_taxonomy(taxid):
    # Function retrieves taxonomy of a hit from NCBI.
    #
    # :param taxid: target Taxonomy ID;
    # :type taxid: int;

    ranks = ('Superkingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

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
            'organism_name': 'NA',
            'Domain': 'NA',
            'Phylum': 'NA',
            'Class': 'NA',
            'Order': 'NA',
            'Family': 'NA',
            'Genus': 'NA',
            'Species': 'NA',
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
    organism_name = re.search(r"\<title\>Taxonomy browser \((.+)\)\</title\>", taxonomy_text).group(1)
    taxonomy_dict['organism_name'] = organism_name

    # Check if species name is specified like other ranks:
    check_direct_species_patt = r'TITLE=\"(species)\"\>([A-Za-z0-9 \.]+)\</a\>'
    match_direct_species = re.search(check_direct_species_patt, taxonomy_text)

    if not match_direct_species is None:
        # If species name is specified like other ranks, merely add it to list:
        taxonomy_dict['Species'] = taxonomy_dict['Genus'] + ' ' + match_direct_species.group(2).partition(' ')[2]
    else:
        # Otherwise we need to parse species name from organism_name

        # Get words
        organism_name_words = organism_name.split(' ')

        # We will take all this words as species name.
        #   Example: MN908947
        try:
            if organism_name_words[1] in second_words_not_species:
                taxonomy_dict['Species'] = taxonomy_dict['Genus'] + ' ' + ' '.join(organism_name_words[1:])
            else:
                taxonomy_dict['Species'] = taxonomy_dict['Genus'] + ' ' + organism_name_words[1]
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

    # Rename 'Superkingdom' -> 'Domain'
    taxonomy_dict['Domain'] = taxonomy_dict['Superkingdom']
    del taxonomy_dict['Superkingdom']

    return taxonomy_dict
# end def


def fill_missing_taxonomy(row):

    if pd.isnull(row['organism_name']):

        print(f'Requesting taxonomy for taxid {row["taxid"]}... ')

        taxonomy_dict = request_missing_taxonomy(row['taxid'])

        # If no taxonomy was retrieved, get it from the RefSeq title
        # Bad way, but no better ways are left
        if taxonomy_dict['organism_name'] == 'NA':
            global asm_sum_df
            seq_title = asm_sum_df[asm_sum_df['asm_acc'] == row['asm_acc']] \
                .reset_index().loc[0, 'title']
            strings_to_rm = (
                ', complete sequence',
                ', complete genome',
                ' map unlocalized',
            )
            for str_to_rm in strings_to_rm:
                seq_title = seq_title.replace(str_to_rm, '')
            # end for
            taxonomy_dict['organism_name'] = seq_title

            print(f'Cannot find taxonomy for taxid {row["taxid"]} at the NCBI website')
            print(f'Using the RefSeq title as the taxonomy name: `{taxonomy_dict["organism_name"]}`')
        else:
            print(taxonomy_dict)
        # end if

        # Fill the taxonomy
        for rank_name, taxon_name in taxonomy_dict.items():
            row[rank_name] = taxon_name
        # end for
    # end if

    return row
# end def



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
    names=[
        'taxid', 'organism_name',
        'Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum',
        'Unknown_1', 'Unknown_2', 'Domain'
    ],
    header=None,
    index_col=False,
    dtype={
        'taxid': pd.Int32Dtype(),
        'organism_name': str,
        'Species': str,
        'Genus': str,
        'Family': str,
        'Order': str,
        'Class': str,
        'Phylum': str,
        'Unknown_1': str,
        'Unknown_2': str,
        'Domain': str
    }
)

# Remove columns of no interest
rankedlineage_df = rankedlineage_df.drop(
    columns=['organism_name', 'Unknown_1', 'Unknown_2'],
    axis=1
)


# Make per-genome taxonomy file

print('Creating taxonomy file')

# Read per-genome taxid file
asm_sum_df = rgIO.read_ass_sum_file(asm_sum_fpath)

# Merge per-genome taxid file to rankedlineage file
taxonomy_df = asm_sum_df.merge(rankedlineage_df, on='taxid', how='left')
del rankedlineage_df

# Sometimes information is missing for some taxids in rankedlineade.dmp
# Request the missing taxonomy from NCBI Taxonomy
missing_taxids = set(
    taxonomy_df[
        pd.isnull(taxonomy_df['organism_name'])
    ]['taxid']
)
if len(missing_taxids) != 0:
    print(f'Taxonomy is missing for {len(missing_taxids)} Taxonomy IDs')
    print('The script will request the taxonomy for them from the NCBI website')
    taxonomy_df = taxonomy_df.apply(fill_missing_taxonomy, axis=1)
# end if
del missing_taxids


# Amend class for Cyanobacteria
taxonomy_df = taxonomy_df.apply(amend_Cyanophyceae, axis=1)
# Amend species names
taxonomy_df = taxonomy_df.apply(fill_empty_species_name, axis=1)

# Order columns
taxonomy_df = taxonomy_df[
    [
        'asm_acc',
        'taxid',
        'organism_name',
        'Species',
        'Genus',
        'Family',
        'Order',
        'Class',
        'Phylum',
        'Domain',
    ]
]


# Write output per-genome file
taxonomy_df.to_csv(
    outfpath,
    sep='\t',
    na_rep='NA',
    header=True,
    index=False,
    encoding='utf-8'
)


print('\nCompleted!')
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
