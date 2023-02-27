#!/usr/bin/env python3

# TODO: add description


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')

import sys
import gzip
import argparse

import numpy as np
import pandas as pd

import src.rg_tools_IO as rgIO


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files


parser.add_argument(
    '-i',
    '--in-ass-sum',
    help="""TSV file (with header) with
    Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '-a',
    '--refseq-catalog',
    help="""A RefSeq "catalog" file of the current release.
This is the file `RefSeq-releaseXXX.catalog.gz` from here:
https://ftp.ncbi.nlm.nih.gov/refseq/release/release-catalog/.
It is better to filter this file with `filter_refseq_catalog.py` before running current script.
""",
    required=True
)

parser.add_argument(
    '-b',
    '--acc-blacklist',
    help="""A TSV file listing RefSeq Accession numbers (without version) to be discarded.
    The file should contain a header.
    Also, it should contains at least one column (of accession numbers).
    The second column (reason for rejection) is optional.""",
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--out-ass-sum',
    help="""a filtered TSV file (with header) with
    GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

args = parser.parse_args()


infpath = os.path.realpath(args.in_ass_sum)
filtered_catalog_fpath = os.path.realpath(args.refseq_catalog)
blacklist_fpath = os.path.realpath(args.acc_blacklist)
outfpath = os.path.realpath(args.out_ass_sum)


# Check existance of the input files

for fpath in (infpath, filtered_catalog_fpath, blacklist_fpath):
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

print(infpath)
print(filtered_catalog_fpath)
print(blacklist_fpath)
print()



ALLOWED_ASM_LEVELS = (
    'Complete Genome',
    'Chromosome'
)


def remove_wgs_assemblies(ass_sum_df):
    return ass_sum_df[
        pd.isnull(ass_sum_df['wgs_master'])
    ].copy()
# end def



def filter_by_ass_level(ass_sum_df):
    global ALLOWED_ASM_LEVELS
    allowed_ass_levels = set(
        map(lambda x: x.upper(), ALLOWED_ASM_LEVELS)
    )
    ass_sum_df['up_assembly_level'] = np.repeat('', ass_sum_df.shape[0])
    ass_sum_df = ass_sum_df.apply(set_up_ass_level, axis=1)

    filt_ass_sum_df = ass_sum_df.query(
        'up_assembly_level in @allowed_ass_levels'
    )

    filt_ass_sum_df = filt_ass_sum_df.drop(columns=['up_assembly_level'])
    return filt_ass_sum_df
# end def

def set_up_ass_level(row):
    row['up_assembly_level'] = row['assembly_level'].upper()
    return row
# end def


def remove_blacklist(ass_sum_df, blacklist_fpath):
    blacklist_accs = read_blacklist(blacklist_fpath)
    ass_sum_df['acc_no_version'] = np.repeat('', ass_sum_df.shape[0])
    ass_sum_df = ass_sum_df.apply(set_acc_no_version, axis=1)

    filt_ass_sum_df = ass_sum_df.query(
        'not acc_no_version in @blacklist_accs'
    )

    filt_ass_sum_df = filt_ass_sum_df.drop(columns=['acc_no_version'])
    return filt_ass_sum_df
# end def

def set_acc_no_version(row):
    row['acc_no_version'] = row['asm_acc'].partition('.')[0]
    return row
# end def

def read_blacklist(blacklist_fpath):
    accessions = set()
    with open(blacklist_fpath, 'rt') as blfile:
        blfile.readline() # pass header
        line = blfile.readline()
        while line != '':
            accessions.add(
                line.split('\t')[0]
            )
            line = blfile.readline()
        # end while
    # end with
    return accessions
# end def


def write_output(ass_sum_df, outfpath):
    with gzip.open(outfpath, 'wt') as outfile:
        ass_sum_df.to_csv(
            outfile,
            sep='\t',
            index=False,
            header=True,
            encoding='utf-8',
            na_rep='NA'
        )
    # end with
# end def



# == Proceed ==

raw_ass_sum_df = rgIO.read_ass_sum_file(infpath, raw_summary=True)
raw_rownum = raw_ass_sum_df.shape[0]
print('Found {:,} genomes at all.'.format(raw_rownum))

# 1.
ass_sum_df = remove_wgs_assemblies(raw_ass_sum_df)
step1_rownum = ass_sum_df.shape[0]
print('1. "Whole genome shotgun" sequences are removed.')
print('   {:,} genomes are retained.'.format(step1_rownum))

# 2.
ass_sum_df = filter_by_ass_level(ass_sum_df)
step2_rownum = ass_sum_df.shape[0]
print('2. Removed assemblies of level other than {}.'.format(ALLOWED_ASM_LEVELS))
print('   {:,} genomes are retained.'.format(step2_rownum))

# 3.
ass_sum_df = remove_blacklist(ass_sum_df, blacklist_fpath)
step3_rownum = ass_sum_df.shape[0]
print('3. Blacklist genomes are removed.')
print('   {:,} genomes are retained.'.format(step3_rownum))

# Output
write_output(ass_sum_df, outfpath)


print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
