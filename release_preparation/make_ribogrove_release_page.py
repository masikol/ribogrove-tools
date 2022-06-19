#!/usr/bin/env python3

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')

import re
import sys
import argparse

import flask
import pandas as pd


from src.ribogrove_size import make_ribogrove_size_dict, format_size_dict
from src.gene_lengths import make_ribogrove_len_dict, format_len_dict
from src.copy_number import make_ribogrove_copy_number_df, format_copy_number_df
from src.top_longest_genes import make_ribogrove_top_longest_df, format_longest_genes_df
from src.top_shortest_genes import make_ribogrove_top_shortest_df, format_shortest_genes_df
from src.top_copy_numbers import make_ribogrove_top_copy_numbers_df, format_top_copy_numbers_df
from src.top_variability import make_ribogrove_top_intragenomic_var_df, format_top_intragenomic_var_df
from src.formatting import format_float_number
from src.strains_names import retrieve_strain_name_en, \
                              retrieve_strain_name_ru, \
                              retrieve_strain_name_ua, \
                              retrieve_strain_name_be

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
    '--raw-fasta',
    help='RiboGrove "raw" fasta file (Archaea + Bacteria)`',
    required=True
)

parser.add_argument(
    '--metadata',
    help='RiboGrove metadata archive (Archaea + Bacteria)`',
    required=True
)

parser.add_argument(
    '--gene-stats-table',
    help='per-gene statistivs for RiboGrove seuences (Archaea + Bacteria)`',
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
raw_fasta_fpath = os.path.abspath(args.raw_fasta)
metadata_fpath = os.path.abspath(args.metadata)
gene_stats_fpath = os.path.abspath(args.gene_stats_table)
entropy_summary_fpath = os.path.abspath(args.entropy_summary)
source_genomes_fpath = os.path.abspath(args.source_genomes)
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
    raw_fasta_fpath,
    metadata_fpath,
    gene_stats_fpath,
    entropy_summary_fpath,
    source_genomes_fpath,
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


# == Some functions ==

def get_file_size_MB(fpath):

    size_in_bytes = os.path.getsize(fpath)
    size_in_megabytes = size_in_bytes / 1024 / 1024

    return size_in_megabytes
# end def get_file_size_MB


# == Proceed ==

# Some Flask stuff
app = flask.Flask('RiboGrove page maker')
app.instance_path = os.path.abspath(os.path.dirname(__file__))

# Get RefSeq release on which the cueernt RiboGrove release is based
refseq_release = ribogrove_release_number.partition('.')[2]

# Calculate file sizes
final_fasta_fsize = get_file_size_MB(final_fasta_fpath)
raw_fasta_fsize = get_file_size_MB(raw_fasta_fpath)
metadata_fsize = get_file_size_MB(metadata_fpath)


# Read per-gene statistics file
gene_stats_df = pd.read_csv(gene_stats_fpath, sep='\t')
# Read info about source genomes
source_genomes_df = pd.read_csv(source_genomes_fpath, sep='\t')
init_columns = gene_stats_df.columns
gene_stats_df = gene_stats_df.merge(
    source_genomes_df[['ass_id', 'title',]],
    on='ass_id',
    how='left'
).drop_duplicates(subset=init_columns)
del init_columns


# Read entropy summary file
entropy_summary_df = pd.read_csv(entropy_summary_fpath, sep='\t')

# RiboGrove size
print('Counting RiboGrove size')
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
print('done\n')

# RiboGrove top shortest genes
print('Finding top shortest RiboGrove genes')
ribogrove_top_shortest_df = make_ribogrove_top_shortest_df(
    gene_stats_df
)
print('done\n')

# RiboGrove top shortest genes
print('Finding top genomes with largest copy numbers RiboGrove genes')
ribogrove_top_copy_numbers_df = make_ribogrove_top_copy_numbers_df(
    gene_stats_df
)
print('done\n')

# RiboGrove top genomes with highest intragenomic variability of target genes
print('Finding top genomes with highest intragenomic variability of target genes')
ribogrove_top_intragenomic_var_df = make_ribogrove_top_intragenomic_var_df(
    entropy_summary_df,
    gene_stats_df
)
print('done\n')


# Language stuff
# (English, Russian)

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
    raw_fasta_fsize_fmt = format_float_number(
        raw_fasta_fsize,
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

    fmt_ribogrove_top_longest_df = format_longest_genes_df(
        ribogrove_top_longest_df,
        thousand_separator,
        decimal_separator
    )

    fmt_ribogrove_top_shortest_df = format_shortest_genes_df(
        ribogrove_top_shortest_df,
        thousand_separator,
        decimal_separator
    )

    fmt_ribogrove_top_copy_numbers_df = format_top_copy_numbers_df(
        ribogrove_top_copy_numbers_df,
        thousand_separator,
        decimal_separator
    )

    fmt_ribogrove_top_intragenomic_var_df = format_top_intragenomic_var_df(
        ribogrove_top_intragenomic_var_df,
        thousand_separator,
        decimal_separator
    )

    # Render the template
    with app.app_context():
        rendered_str = flask.render_template(
                template_fpath,
                ribogrove_release_number=ribogrove_release_number,
                refseq_release=refseq_release,
                ribogrove_release_date=ribogrove_release_date,
                final_fasta_fsize_fmt=final_fasta_fsize_fmt,
                raw_fasta_fsize_fmt=raw_fasta_fsize_fmt,
                metadata_fsize_fmt=metadata_fsize_fmt,
                ribogrove_size_dict=fmt_ribogrove_size_dict,
                ribogrove_len_dict=fmt_ribogrove_len_dict,
                ribogrove_copy_number_df=fmt_ribogrove_copy_number_df,
                ribogrove_top_longest_df=fmt_ribogrove_top_longest_df,
                ribogrove_top_shortest_df=fmt_ribogrove_top_shortest_df,
                ribogrove_top_copy_numbers_df=fmt_ribogrove_top_copy_numbers_df,
                ribogrove_top_intragenomic_var_df=fmt_ribogrove_top_intragenomic_var_df,
                retrieve_strain_name=retrieve_strain_name
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
