#!/usr/bin/env python3

# The script merges two files together:
# 1) a TSV file, which is output of the script `assembly2gi_numbers.py` (`-c` option);
# 2) a TSV file, which is output of the script `gis_to_accs.py` (`-c` option).

# The output file (`-o` option) is a TSV file of 4 columns (`ass_id`, `gi_number`, `acc`, `title`).
#   In this file, every line corresponds to a single RefSeq record, and the line
#   contains Assembly ID (`ass_id`), RefSeq GI number (`gi_number`), ACCESSION.VERSION (`acc`),
#   and the title of the record (`title`).

## Command line arguments
### Input files:
# 1. `-s / --assm-2-gi-file` -- input TSV file mapping Assembly IDs to RefSeq GI numbers. Mandatory.
# 2. `-c / --gi-2-acc-file` -- input TSV file mapping RefSeq GI numbers to RefSeq ACCESSION.VERSIONs
#   and titles. Mandatory.

### Output files:
# 1. `-o / --outfile` -- output TSV files where Assembly IDs are mapped to ACCESSION.VERSIONs
#   and RefSeq Titles. Mandatory.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import gzip
import argparse


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--raw-refseq-catalog',
    help="""A RefSeq "catalog" file of the current release.
This is the file `RefSeq-releaseXXX.catalog.gz` from here:
https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/""",
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='filtered RefSeq catalog: only genomic sequences, only bacteria and archaea',
    required=True
)

args = parser.parse_args()

raw_refseq_catalog_fpath = os.path.realpath(args.raw_refseq_catalog)
outfpath = os.path.realpath(args.outfile)


# Check existance of the input files

for fpath in (raw_refseq_catalog_fpath,):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# end for

if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

print(raw_refseq_catalog_fpath)
print()


# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
non_genomic_prefixes = {
    'NM_', 'NR_', 'XM_', 'XR_', 'AP_', 'NP_', 'YP_', 'XP_', 'WP_',
}

target_organisms = ('bacteria', 'archaea')

def is_target_organism(dir_column_val):
    dir_column_val = dir_column_val.lower()

    return any(
        map(
            lambda o: o in dir_column_val,
            target_organisms
        )
    )
# end def is_target_organism


# == Proceed ==


# Read the catalog file
if raw_refseq_catalog_fpath.endswith('.gz'):
    open_func = gzip.open
else:
    open_func = open
# end if

print(f'Filtering the catalog file `{raw_refseq_catalog_fpath}`...')


status_increment = 100000
next_done_num = status_increment

acc_column_index = 2
dir_column_index = 3
separator = '\t'

kept_number = 0

with open_func(raw_refseq_catalog_fpath, 'rt') as catalog_file, \
     gzip.open(outfpath, 'wt') as filtered_catalog_file:

    print('0 lines processed', end=' '*10)

    for i, line in enumerate(catalog_file):

        if i == next_done_num:
            print(
                '\r{:,} lines processed; {:,} lines are passing' \
                    .format(i, kept_number),
                end=' '*10
            )
            next_done_num += status_increment
        # end if

        line_vals = line.split(separator)
        prefix = line_vals[acc_column_index][:3]

        if prefix in non_genomic_prefixes:
            continue
        # end if
        if not is_target_organism(line_vals[dir_column_index]):
            continue
        # end if

        kept_number += 1
        filtered_catalog_file.write(line)
    # end for
# end with

print(
    '\r{:,} lines processed; {:,} lines passed' \
        .format(i+1, kept_number)
)
print('done\n')


print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
