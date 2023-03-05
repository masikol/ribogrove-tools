#!/usr/bin/env python3

# The script assigns categories to downloaded genomes. Categories are assigned according
#   to the reliability of a genome assembly.
# The categories are the following:
# Category 1. A genome is not of category 3, and it was sequenced using PacBio or ONT+Illumina.
# Category 2. A genome is not of category 3, and it was sequenced neither using PacBio nor ONT+Illumina.
# Category 3. At least one of the following is true:
#   - A genome has at least one degenerate base in its SSU gene sequences.
#   - At least one of the genomic sequences contains phrase "map unlocalized" in it title,
#     and the sequence contains an SSU gene (or a part of it).

## Command line arguments

### Input files:
# 1. `-f / --all-fasta-file` -- a fasta file with all extracted genes sequences.
#   This file is the output of the script `extract_16S.py`.
#   Mandatory.
# 2. `-s / --all-stats-file` -- a file with per-replicon statistics of input 16S gene sequences.
#   This file is the output of the script `extract_16S.py`, too.
#   Mandatory.
# 3. `-g / --genomes-dir` -- the directory where the downloaded genome files are located
#   (see script `download_genomes.py`).
#   Mandatory.

### Output files:
# 1. `-o / --outfile` -- output file mapping genes Assembly accession numbers to categories.
#   Mandatory.

### Dependencies:
# 1. `--seqkit` -- a `seqkit` executable: github.com/shenwei356/seqkit. Mandatory.


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
    '-f',
    '--all-fasta-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-s',
    '--all-stats-file',
    help="""TSV file (with header) containing per-replicons SSU gene statistics
    reported by extract_16S.py""",
    required=True
)

parser.add_argument(
    '-g',
    '--genomes-dir',
    help='directory that contains downloaded gbk.gz files',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='output file mapping genes Assembly accession numbers to categories',
    required=True
)

# Dependencies

parser.add_argument(
    '--seqkit',
    help='seqkit executable',
    required=True
)


args = parser.parse_args()


# == Import them now ==
import re
import sys
import gzip
import subprocess as sp
from typing import Dict, List, Sequence, TextIO

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from src.ribogrove_seqID import update_seqID_column_v2_to_v3
from src.file_navigation import get_asm_report_fpath


# For convenience
fasta_seqs_fpath = os.path.abspath(args.all_fasta_file)
in_stats_fpath = os.path.abspath(args.all_stats_file)
genomes_dirpath = os.path.abspath(args.genomes_dir)
outfpath = os.path.abspath(args.outfile)
seqkit_fpath = os.path.abspath(args.seqkit)


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, in_stats_fpath, seqkit_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Check if input directory exists
if not os.path.exists(genomes_dirpath):
    print(f'Error: directory `{genomes_dirpath}` does not exist!')
    sys.exit(1)
# end if

# Check if seqkit executable is actually executable
if not os.access(seqkit_fpath, os.X_OK):
    print(f'Error: file `{seqkit_fpath}` is not executable!')
    sys.exit(1)
# end if

# Create output directories if needed
if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

print(fasta_seqs_fpath)
print(in_stats_fpath)
print(genomes_dirpath)
print(seqkit_fpath)
print()


# Paths to files containing marker words for identifying sequencing technologies
data_dirpath = os.path.join(
    os.path.dirname(__file__),
    'data'
)
pacbio_vocab_fpath = os.path.realpath(
    os.path.join(data_dirpath, 'seqtech_dicts', 'pacbio')
)
illumina_vocab_fpath = os.path.realpath(
    os.path.join(data_dirpath, 'seqtech_dicts', 'illumina')
)
nanopore_vocab_fpath = os.path.realpath(
    os.path.join(data_dirpath, 'seqtech_dicts', 'ont')
)
del data_dirpath


# These keywords might be recognized as "ONT", but they are not
#   releated to Oxford Nanopore.
# These words must be removed from seqtech string befor searching for keywords
seem_like_ont_but_not = {
    'IONTORRENT',
    # ~~~
    'CONTIG', # see NC_020549.1
    # ~~~
}

# Let them be str from start, not int
CATEGORY_1 = '1'
CATEGORY_2 = '2'
CATEGORY_3 = '3'

SEQTECH_PATTERN = re.compile(
    r'Sequencing technology:([^#]+)'
)



def find_degenerate_in_16S(fasta_seqs_fpath: str, seqkit_fpath: str) -> List[str]:
    # Function reports Assembly IDs of genomes which contain degenerate bases in their SSU genes.

    # Configure command reporting accessions of sequences which contain degenerate bases in their SSU genes
    cmd = ' '.join(
        [
            '{} grep -srp "[RYWSKMHVBDN]" {}'.format(seqkit_fpath, fasta_seqs_fpath),
            ' | {} seq -ni'.format(seqkit_fpath),
            ' | cut -f1 -d":"',
            ' | sort',
            ' | uniq'
        ]
    )
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout, stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error at popen finding degenerate in 16S genes')
        print(stderr.decode('utf-8'))
        sys.exit(1)
    else:
        # Parse accessions
        asm_accs_degen_in_16S = set(stdout.decode('utf-8').split('\n'))
    # end if

    return asm_accs_degen_in_16S
# end def


def read_seqtech_vocab(fpath: str) -> Sequence[str]:
    # Function makes vocabulary of marker seqtech words
    with open(fpath, 'rt') as vocab_file:
        vocab = tuple(
            map(
                lambda s: s.strip().upper(),
                vocab_file.readlines()
            )
        )
    # end with
    return vocab
# end def


def is_pacbio(seqtech_str: str) -> bool:
    # Function decides if genome was sequenced with PacBio

    global pacbio_vocab

    return any(
        tuple(
            map(
                lambda x: x in seqtech_str,
                pacbio_vocab
            )
        )
    )
# end def

def is_illumina(seqtech_str: str) -> bool:
    # Function decides if genome was sequenced with Illumina

    global illumina_vocab

    return any(
        tuple(
            map(
                lambda x: x in seqtech_str,
                illumina_vocab
            )
        )
    )
# end def

def is_nanopore(seqtech_str: str) -> bool:
    # Function decides if genome was sequenced with Oxford Nanopore

    global nanopore_vocab
    global seem_like_ont_but_not

    for keyword in nanopore_vocab:
        if keyword != 'ONT':
            # If keyword is not "ONT" we proceed just like with any other keyword
            if keyword in seqtech_str:
                return True
            # end if
        else:
            # It kwyword is "ONT", we should previously remove words from `seem_like_ont_but_not`
            #   before searching the keyword
            for word in seem_like_ont_but_not:
                seqtech_str = seqtech_str.replace(word, '')
            # end for
            if keyword in seqtech_str:
                return True
            # end if
        # end if
    # end for
    return False
# end def


def parse_seqtech(asm_acc, genomes_dirpath):
    global SEQTECH_PATTERN

    asm_report_fpath = get_asm_report_fpath(asm_acc, genomes_dirpath)
    with open(asm_report_fpath, 'rt') as report_file:
        seqtech_re_obj = re.search(SEQTECH_PATTERN, report_file.read())
    # end with

    if seqtech_re_obj is None:
        return None
    else:
        return seqtech_re_obj.group(1).strip().upper()
    # end if
# end def



# Read per-replicon statistics
stats_df = pd.read_csv(
    in_stats_fpath,
    sep='\t'
)


# Make vocabularies of seqtech keywords
pacbio_vocab   = read_seqtech_vocab(pacbio_vocab_fpath)
illumina_vocab = read_seqtech_vocab(illumina_vocab_fpath)
nanopore_vocab = read_seqtech_vocab(nanopore_vocab_fpath)


# Get Assembly IDs of genomes having degenerate bases in SSU genes
print('Searching for 16S genes containing degenerate bases...')
asm_accs_degen_in_16S = find_degenerate_in_16S(fasta_seqs_fpath, seqkit_fpath)
print(
    'Found {:,} genomes containing 16S genes with degenerate bases\n' \
        .format(len(asm_accs_degen_in_16S))
)


# == Proceed ==

print('Starting assigning categories for genomes')

with open(outfpath, 'wt') as outfile:

    # Write header
    outfile.write('asm_acc\tcategory\tseqtech\tdegenerate_in_16S\tunlocalized_16S\n')

    # Get all Assembly accessions
    all_asm_accs = tuple(set(stats_df['asm_acc']))

    status_step = 50
    next_status = min(status_step, len(all_asm_accs))

    sys.stdout.write('0/{} genomes are done'.format(len(all_asm_accs)))
    sys.stdout.flush()

    # Iterate over Assembly IDs
    for i, asm_acc in enumerate(all_asm_accs):
        # Genome has (maybe, patrial) SSU genes in "map unlocalized" sequences
        unlocalized_16S = False
        # Genome has degenerate bases in SSU genes
        degenerate_in_16S = asm_acc in asm_accs_degen_in_16S

        # Sequencing technology
        seqtech = parse_seqtech(asm_acc, genomes_dirpath)

        # Get rows corresponding to current assembly
        curr_asm_df = stats_df[stats_df['asm_acc'] == asm_acc]

        # Iterate over ACCESSION.VERSION's of current genome
        for _, row in curr_asm_df.iterrows():
            asm_acc = row['asm_acc']
            title   = row['title']

            # Update flag `unlocalized_16S`
            map_unlocalized = 'MAP UNLOCALIZED' in title.upper()
            unlocalized_16S = unlocalized_16S or (map_unlocalized and row['num_genes'] != 0)
        # end for

        # Assign category to the genome
        category = None
        if degenerate_in_16S or unlocalized_16S:
            category = CATEGORY_3
        elif not seqtech is None:
            seqtech_is_cat1 = is_pacbio(seqtech) \
                              or (is_illumina(seqtech) and is_nanopore(seqtech))
            if seqtech_is_cat1:
                category = CATEGORY_1
            else:
                category = CATEGORY_2
            # end if
        else:
            category = CATEGORY_2
        # end if

        # Write ouput
        out_row_str = '\t'.join(
            [
                asm_acc, category,
                seqtech if not seqtech is None else 'NA',
                '1' if degenerate_in_16S else '0',
                '1' if   unlocalized_16S else '0',
            ]
        )
        outfile.write('{}\n'.format(out_row_str))

        if (i+1) == next_status:
            sys.stdout.write(
                '\r{}/{} genomes are done'.format(i+1, len(all_asm_accs))
            )
            sys.stdout.flush()
            next_status = min(
                next_status + status_step,
                len(all_asm_accs)
            )
        # end if
    # end for
# end with

print('\n\nCompleted!')
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
