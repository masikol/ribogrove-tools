#!/usr/bin/env python3

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')

import re
import sys
import json
import argparse
from functools import reduce
from collections import OrderedDict

import flask
import numpy as np
import pandas as pd


from src.ribogrove_size import make_ribogrove_size_dict, format_size_dict
from src.gene_lengths import make_ribogrove_len_dict, format_len_dict
from src.copy_number import make_ribogrove_copy_number_df, format_copy_number_df
from src.top_longest_genes import make_ribogrove_top_longest_df, format_longest_genes_df
from src.top_shortest_genes import make_ribogrove_top_shortest_df, format_shortest_genes_df
from src.top_copy_numbers import make_ribogrove_top_copy_numbers_df, format_top_copy_numbers_df
from src.top_variability import make_ribogrove_top_intragenomic_var_df, format_top_intragenomic_var_df
from src.primers_coverage import make_ribogrove_primer_coverage_df, format_primer_coverage_df
from src.formatting import format_float_number
from src.strains_names import retrieve_strain_name_en, \
                              retrieve_strain_name_ru, \
                              retrieve_strain_name_ua, \
                              retrieve_strain_name_be, \
                              italicize_candidatus

# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input data

parser.add_argument(
    '-r',
    '--release-num',
    help='RiboGrove release number, e.g. `2.208`',
    required=True
)

parser.add_argument(
    '-d',
    '--release-date',
    help='RiboGrove release date, e.g. `2021-10-28`',
    required=True
)

parser.add_argument(
    '--final-fasta',
    help='RiboGrove final fasta file (Archaea + Bacteria)`',
    required=True
)

parser.add_argument(
    '--metadata',
    help='RiboGrove metadata archive (Archaea + Bacteria)`',
    required=True
)

parser.add_argument(
    '--base-counts',
    help='per-gene base counts table for RiboGrove seuences (Archaea + Bacteria)`',
    required=True
)

parser.add_argument(
    '--taxonomy',
    help='taxonomy table for RiboGrove seuences (Archaea + Bacteria)`',
    required=True
)

parser.add_argument(
    '--categories',
    help='categories table for RiboGrove seuences (Archaea + Bacteria)`',
    required=True
)

parser.add_argument(
    '--entropy-summary',
    help='per-genome summary of intragenomic variablility (of entropy) (Archaea + Bacteria)`',
    required=True
)

parser.add_argument(
    '--source-genomes',
    help='a file of information about what genomes were used for the RiboGrove construction',
    required=True
)

parser.add_argument(
    '--primers-cov',
    help='the file of primer coverage summary, `primer_pair_genomic_coverage.tsv`',
    required=True
)

parser.add_argument(
    '--archive',
    help='make pages for an archive realease, i.e. they shall be reduced',
    action='store_true'
)

# Output files

parser.add_argument(
    '-o',
    '--outdir',
    help='output directory for rendered HTML files in all available languages',
    required=True
)

# Dependencies

parser.add_argument(
    '--seqkit',
    help='path to seqkit executable`',
    required=True
)


args = parser.parse_args()

ribogrove_release_number = args.release_num
ribogrove_release_date = args.release_date
final_fasta_fpath = os.path.abspath(args.final_fasta)
metadata_fpath = os.path.abspath(args.metadata)
base_counts_fpath = os.path.abspath(args.base_counts)
taxonomy_fpath = os.path.abspath(args.taxonomy)
categories_fpath = os.path.abspath(args.categories)
entropy_summary_fpath = os.path.abspath(args.entropy_summary)
source_genomes_fpath = os.path.abspath(args.source_genomes)
primers_fpath = os.path.abspath(args.primers_cov)
archive = args.archive
outdpath = os.path.abspath(args.outdir)
seqkit_fpath = os.path.abspath(args.seqkit)



# == Validate the arguments ==

if re.match(r'^[\.0-9]+$', ribogrove_release_number) is None:
    print('Invalid format of RiboGrove release number!')
    print(f'Must be like `2.208`. You entered `{ribogrove_release_number}`')
    sys.exit(1)
# end if

if re.match(r'^[\-0-9]+$', ribogrove_release_date) is None:
    print('Invalid format of RiboGrove release date!')
    print(f'Must be like `2021-10-28`. You entered `{ribogrove_release_date}`')
    sys.exit(1)
# end if

input_fpaths = (
    final_fasta_fpath,
    metadata_fpath,
    base_counts_fpath,
    taxonomy_fpath,
    categories_fpath,
    entropy_summary_fpath,
    source_genomes_fpath,
    primers_fpath
)

for fpath in input_fpaths:
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist')
        sys.exit(1)
    # end if
# end for

if not os.path.isdir(outdpath):
    try:
        os.makedirs(outdpath)
    except OSError as err:
        print(f'Error: cannot create directory `{outdpath}`')
        print(str(err))
        sys.exit(1)
    # end try
# end if

del input_fpaths, fpath


STRAIN_DESIGNATION_PATTERN = re.compile(
    r'strain=(.+)'
)


# == Some functions ==

def make_gene_stats_df(base_counts_fpath,
                       taxonomy_fpath,
                       categories_fpath):
    base_counts_df = pd.read_csv(base_counts_fpath, sep='\t')
    taxonomy_df = pd.read_csv(taxonomy_fpath, sep='\t')
    categories_df = pd.read_csv(categories_fpath, sep='\t')

    base_counts_df['asm_acc'] = np.repeat('', base_counts_df.shape[0])
    base_counts_df = base_counts_df.apply(set_asm_acc, axis=1)

    gene_stats_df = base_counts_df[['asm_acc', 'seqID', 'len']].merge(
        taxonomy_df,
        on='asm_acc',
        how='left'
    ).merge(
        categories_df[['asm_acc', 'category']],
        on='asm_acc',
        how='left'
    )
    return gene_stats_df
# end def

def set_asm_acc(row):
    # TODO: parse asm_acc using the function in src...
    row['asm_acc'] = row['seqID'].partition(':')[0]
    return row
# end def

def get_file_size_MB(fpath):

    size_in_bytes = os.path.getsize(fpath)
    size_in_megabytes = size_in_bytes / 1024 / 1024

    return size_in_megabytes
# end def


def set_strain_name(row):
    global STRAIN_DESIGNATION_PATTERN
    if pd.isnull(row['infraspecific_name']):
        row['infraspecific_name'] = ''
    # end if
    strain_name_reobj = STRAIN_DESIGNATION_PATTERN.search(row['infraspecific_name'])
    if not strain_name_reobj is None:
        strain_designation = strain_name_reobj.group(1)
        if not row['organism_name'].endswith(strain_designation):
            row['strain_name'] = '{} strain {}'.format(
                row['organism_name'],
                strain_designation
            )
        else:
            row['strain_name'] = row['organism_name']
        # end if
    else:
        row['strain_name'] = row['organism_name']
    # end if
    return row
# end def

def parse_primer_pairs():
    # TODO: deduplicate code
    primers_pairs_fpath = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'create_RiboGrove', 'collect_and_filter', 'scripts', 'data', 'primers', 'primer_pairs.json'
    )
    with open(primers_pairs_fpath, 'rt') as infile:
        primer_pairs = json.load(infile)
    # end with

    bacterial_primer_pairs = transform_primer_pair_dict(primer_pairs['Bacteria'])
    archaeal_primer_pairs  = transform_primer_pair_dict(primer_pairs['Archaea'])

    unwanted_primer_pairs = parse_unwanted_primer_pairs()

    for pair_name in unwanted_primer_pairs:
        if pair_name in bacterial_primer_pairs:
            del bacterial_primer_pairs[pair_name]
        # end if
        if pair_name in archaeal_primer_pairs:
            del archaeal_primer_pairs[pair_name]
        # end if
    # end for

    return bacterial_primer_pairs, archaeal_primer_pairs
# end def

def transform_primer_pair_dict(primer_pairs):
    all_primer_pair_dict = OrderedDict()
    for nameF, nameR, v_region_name in primer_pairs:
        all_primer_pair_dict['{}-{}'.format(nameF, nameR)] = v_region_name
    # end for
    return all_primer_pair_dict
# end def

def parse_unwanted_primer_pairs():
    unwanted_primer_pairs_fpath = os.path.join(
        os.path.dirname(os.path.abspath(__file__)),
        'unwanted_primer_pairs.txt'
    )
    with open(unwanted_primer_pairs_fpath, 'rt') as input_handle:
        primer_pairs = tuple(
            map(
                str.strip,
                input_handle.readlines()
            )
        )
    # end with
    return primer_pairs
# end def


def split_top_df(df, top_n=10):
    bacteria_df = df.query('Domain == "Bacteria"')
    archaea_df = df.query('Domain == "Archaea"')

    really_top_df = pd.concat(
        [
            bacteria_df.iloc[:top_n,],
            archaea_df.iloc[:top_n,],
        ]
    )
    rest_df = pd.concat(
        [
            bacteria_df.iloc[top_n:,],
            archaea_df.iloc[top_n:,],
        ]
    )

    return really_top_df, rest_df
# end def


# == Proceed ==

# Some Flask stuff
app = flask.Flask('RiboGrove page maker')
app.instance_path = os.path.abspath(os.path.dirname(__file__))

# Get RefSeq release on which the cueernt RiboGrove release is based
refseq_release = ribogrove_release_number.partition('.')[2]

# Calculate file sizes
final_fasta_fsize = get_file_size_MB(final_fasta_fpath)
metadata_fsize = get_file_size_MB(metadata_fpath)


# Read per-gene statistics file
gene_stats_df = make_gene_stats_df(
    base_counts_fpath,
    taxonomy_fpath,
    categories_fpath
)
# Read info about source genomes
source_genomes_df = pd.read_csv(source_genomes_fpath, sep='\t')
source_genomes_df['strain_name'] = np.repeat('', source_genomes_df.shape[0])
source_genomes_df = source_genomes_df.apply(set_strain_name, axis=1)

# Combine the two
init_columns = gene_stats_df.columns
gene_stats_df = gene_stats_df.merge(
    source_genomes_df[['asm_acc', 'strain_name',]],
    on='asm_acc',
    how='left'
).drop_duplicates(subset=init_columns)
del init_columns

# Parse primer pair data
bacterial_primer_pairs, archaeal_primer_pairs = parse_primer_pairs()

# Read entropy summary file
entropy_summary_df = pd.read_csv(entropy_summary_fpath, sep='\t')

# RiboGrove size
print('Counting RiboGrove sequences')
ribogrove_size_dict = make_ribogrove_size_dict(
    final_fasta_fpath,
    gene_stats_df,
    seqkit_fpath
)
print('done\n')

# RiboGrove gene lengths
print('Counting RiboGrove gene lengths')
ribogrove_len_dict = make_ribogrove_len_dict(
    gene_stats_df
)
print('done\n')

# RiboGrove copy number
print('Counting RiboGrove copy numbers')
ribogrove_copy_number_df = make_ribogrove_copy_number_df(
    gene_stats_df
)
print('done\n')

# RiboGrove top longest genes
print('Finding top longest RiboGrove genes')
ribogrove_top_longest_df = make_ribogrove_top_longest_df(
    gene_stats_df
)
really_top_longest_df, rest_top_longest_df = split_top_df(ribogrove_top_longest_df)
print('done\n')

# RiboGrove top shortest genes
print('Finding top shortest RiboGrove genes')
ribogrove_top_shortest_df = make_ribogrove_top_shortest_df(
    gene_stats_df
)
really_top_shortest_df, rest_top_shortest_df = split_top_df(ribogrove_top_shortest_df)
print('done\n')

# RiboGrove top shortest genes
print('Finding top genomes with largest copy numbers RiboGrove genes')
ribogrove_top_copy_numbers_df = make_ribogrove_top_copy_numbers_df(
    gene_stats_df
)
really_top_gcn_df, rest_top_gcn_df = split_top_df(ribogrove_top_copy_numbers_df)
print('done\n')

# RiboGrove top genomes with highest intragenomic variability of target genes
print('Finding top genomes with highest intragenomic variability of target genes')
ribogrove_top_intragenomic_var_df = make_ribogrove_top_intragenomic_var_df(
    entropy_summary_df,
    gene_stats_df
)
print('done\n')

# PCR primer coverage per phylum
print('Calculating PCR primer coverage per phylum')
ribogrove_primer_coverage_df = make_ribogrove_primer_coverage_df(
    primers_fpath
)
print('done\n')


# Language stuff
# (English, Belarusian, Ukrainian, Russian)

template_fpaths = (
    'ribogrove_page_template_en.tpl',
    'ribogrove_page_template_ru.tpl',
    'ribogrove_page_template_ua.tpl',
    'ribogrove_page_template_be.tpl',
)

thousand_separators = (
    ',',
    '&nbsp;',
    '&nbsp;',
    '&nbsp;',
)

decimal_separators = (
    '.',
    ',',
    ',',
    ',',
)

strains_names_functions = (
    retrieve_strain_name_en,
    retrieve_strain_name_ru,
    retrieve_strain_name_ua,
    retrieve_strain_name_be,
)

# Paths to output files
rendered_html_fpaths = tuple(
    map(
        lambda f: os.path.join(
            outdpath,
            f.replace('template', ribogrove_release_number).replace('.tpl', '.html')
        ),
        template_fpaths
    )
)

language_zip = zip(
    template_fpaths,
    thousand_separators,
    decimal_separators,
    rendered_html_fpaths,
    strains_names_functions
)


# Render the templates for all languages

print('\nRendering HTML pages:')

for template_fpath, thousand_separator, decimal_separator, outfpath, retrieve_strain_name in language_zip:

    # Format numeric data

    final_fasta_fsize_fmt = format_float_number(
        final_fasta_fsize,
        thousand_separator,
        decimal_separator,
        2
    )
    metadata_fsize_fmt = format_float_number(
        metadata_fsize,
        thousand_separator,
        decimal_separator,
        2
    )

    fmt_ribogrove_size_dict = format_size_dict(
        ribogrove_size_dict,
        thousand_separator,
        decimal_separator
    )

    fmt_ribogrove_len_dict = format_len_dict(
        ribogrove_len_dict,
        thousand_separator,
        decimal_separator
    )

    fmt_ribogrove_copy_number_df = format_copy_number_df(
        ribogrove_copy_number_df,
        thousand_separator,
        decimal_separator
    )

    fmt_really_top_longest_df = format_longest_genes_df(
        really_top_longest_df,
        thousand_separator,
        decimal_separator
    )
    fmt_rest_top_longest_df = format_longest_genes_df(
        rest_top_longest_df,
        thousand_separator,
        decimal_separator
    )

    fmt_really_top_shortest_df = format_shortest_genes_df(
        really_top_shortest_df,
        thousand_separator,
        decimal_separator
    )
    fmt_rest_top_shortest_df = format_shortest_genes_df(
        rest_top_shortest_df,
        thousand_separator,
        decimal_separator
    )

    fmt_really_top_gcn_df = format_top_copy_numbers_df(
        really_top_gcn_df,
        thousand_separator,
        decimal_separator
    )
    fmt_rest_top_gcn_df = format_top_copy_numbers_df(
        rest_top_gcn_df,
        thousand_separator,
        decimal_separator
    )

    fmt_ribogrove_top_intragenomic_var_df = format_top_intragenomic_var_df(
        ribogrove_top_intragenomic_var_df,
        thousand_separator,
        decimal_separator
    )

    fmt_ribogrove_primers_cov_df = format_primer_coverage_df(
        ribogrove_primer_coverage_df,
        thousand_separator,
        decimal_separator
    )

    # Render the template
    with app.app_context():
        rendered_str = flask.render_template(
                template_fpath,
                archive=archive,
                ribogrove_release_number=ribogrove_release_number,
                refseq_release=refseq_release,
                ribogrove_release_date=ribogrove_release_date,
                final_fasta_fsize_fmt=final_fasta_fsize_fmt,
                metadata_fsize_fmt=metadata_fsize_fmt,
                ribogrove_size_dict=fmt_ribogrove_size_dict,
                ribogrove_len_dict=fmt_ribogrove_len_dict,
                ribogrove_copy_number_df=fmt_ribogrove_copy_number_df,
                really_top_longest_df=fmt_really_top_longest_df,
                rest_top_longest_df=fmt_rest_top_longest_df,
                really_top_shortest_df=fmt_really_top_shortest_df,
                rest_top_shortest_df=fmt_rest_top_shortest_df,
                really_top_gcn_df=fmt_really_top_gcn_df,
                rest_top_gcn_df=fmt_rest_top_gcn_df,
                ribogrove_top_intragenomic_var_df=fmt_ribogrove_top_intragenomic_var_df,
                retrieve_strain_name=retrieve_strain_name,
                ribogrove_primers_cov_df=fmt_ribogrove_primers_cov_df,
                italicize_candidatus=italicize_candidatus,
                bacterial_primer_pairs=bacterial_primer_pairs,
                archaeal_primer_pairs=archaeal_primer_pairs,
                unwanted_primer_pairs=', '.join(parse_unwanted_primer_pairs())
        )
    # end with

    rendered_str = rendered_str.replace('\n\n', '\n')

    # Write rendered HTML to the output file
    with open(outfpath, 'wt') as outfile:
        outfile.write(rendered_str)
    # end with
    print(outfpath)
# end for

print('Completed!')
print(outdpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
