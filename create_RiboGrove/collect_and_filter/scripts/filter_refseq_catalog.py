#!/usr/bin/env python3

# The script filters a RefSeq .catalog.gz file
#   (specifically, file `RefSeq-releaseXXX.catalog.gz`
#   from `https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/`):

# The script reads from stdin and writes to stdout.


import os
import sys
from src.rg_tools_time import get_time

sys.stderr.write(
    '\n|=== {} STARTING SCRIPT `{}` ===|\n\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)


# https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
NON_GENOMIC_PREFIXES = {
    'NM_', 'NR_', 'XM_', 'XR_', 'AP_', 'NP_', 'YP_', 'XP_', 'WP_',
}

TARGET_ORGANISMS = ('bacteria', 'archaea')

def is_target_organism(dir_column_val):
    dir_column_val = dir_column_val.lower()

    return any(
        map(
            lambda o: o in dir_column_val,
            TARGET_ORGANISMS
        )
    )
# end def


# == Proceed ==

ACC_COLUMN_INDEX = 2
DIR_COLUMN_INDEX = 3
SEPARATOR = '\t'

kept_number = 0

for i, line in enumerate(sys.stdin):

    line_vals = line.split(SEPARATOR)
    prefix = line_vals[ACC_COLUMN_INDEX][:3]

    if prefix in NON_GENOMIC_PREFIXES:
        continue
    # end if
    if not is_target_organism(line_vals[DIR_COLUMN_INDEX]):
        continue
    # end if

    kept_number += 1
    sys.stdout.write(line)
# end for

sys.stderr.write(
    '\r{:,} lines processed; {:,} lines passed\n' \
        .format(i+1, kept_number)
)
sys.stderr.write('done\n\n')


sys.stderr.write('\nCompleted!\n')
sys.stderr.write(
    '\n\n|=== {} EXITTING SCRIPT `{}` ===|\n\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
