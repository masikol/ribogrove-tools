#!/usr/bin/env python3

# The script maps Assembly IDs to TaxIDs using elink utility (https://www.ncbi.nlm.nih.gov/books/NBK25497/)

# Input files:
# 1. -i/--assm-acc-file -- is output of script merge_assID2acc_and_remove_WGS.py.
#   It has 4 columns: ass_id, gi_number, acc, title.
# 2. -f/--all-fasta-file -- fasta file of SSU gene sequences

# Output files:
# 1. --per-genome-outfile -- output file mapping Assembly IDs to taxIDs

# Dependencies:
# 1. --seqkit seqkit executable

print(f'\n|=== STARTING SCRIPT `{__file__}` ===|\n')


import os
import sys
import time
import argparse
import subprocess as sp
from typing import List, Dict

import pandas as pd
from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'


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
    '-f',
    '--all-fasta-file',
    help='input fasta file of SSU gene sequences',
    required=True
)

# Output files

parser.add_argument(
    '--per-genome-outfile',
    help='output file mapping Assembly IDs to taxIDs',
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
assm_acc_fpath = os.path.abspath(args.assm_acc_file)
fasta_seqs_fpath = os.path.abspath(args.all_fasta_file)
per_genome_outfpath = os.path.abspath(args.per_genome_outfile)
seqkit_fpath = os.path.abspath(args.seqkit)


# Check existance of all input files and dependencies
for fpath in (assm_acc_fpath, fasta_seqs_fpath, seqkit_fpath):
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
if not os.path.isdir(os.path.dirname(per_genome_outfpath)):
    try:
        os.makedirs(os.path.dirname(per_genome_outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(per_genome_outfpath)}`')
        sys.exit(1)
    # end try
# end if


print(assm_acc_fpath)
print(fasta_seqs_fpath)
print(seqkit_fpath)
print()


def get_genes_seqIDs(fasta_seqs_fpath: str, seqkit_fpath: str) -> List[str]:
    # Function reports all seqIDs of sequences from given fasta file fasta_seqs_fpath.

    # Configure command reporting seqIDs of fasta file
    cmd = f'{seqkit_fpath} seq -ni {fasta_seqs_fpath}'
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error at popen: extracting genes\' seqIDs')
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(pipe.returncode)
    else:
        # Parse seqIDs
        genes_seqIDs = list(stdout_stderr[0].decode('utf-8').split('\n'))
    # end if

    return genes_seqIDs
# end def get_genes_seqIDs


def make_acc_seqIDs_dict(fasta_seqs_fpath: str, seqkit_fpath: str) -> Dict[str, List[str]]:
    # Function creates dictionary that maps accessions to seqIDs

    # Get all seqIDs of gene sequences.
    # We want to have our seqIDs in the same order as they exist in original fasta file,
    #   therefore, we do this `reversed` here
    genes_seqIDs = list(
        reversed(
            get_genes_seqIDs(fasta_seqs_fpath, seqkit_fpath)
        )
    )

    acc_seqIDs_dict = dict()

    for _ in range(len(genes_seqIDs)):

        seqID = genes_seqIDs.pop() # get next seqID
        acc = seqID.partition(':')[0] # parse ACCESSION.VERSION from seqID

        # Fill dictionary
        try:
            acc_seqIDs_dict[acc].append(seqID)
        except KeyError:
            acc_seqIDs_dict[acc] = [seqID]
        # end try
    # end for

    return acc_seqIDs_dict
# end def make_acc_seqIDs_dict


# Read Assembly IDs, ACCESSION.VERSION's dataframe
assm_acc_df = pd.read_csv(
    assm_acc_fpath,
    sep='\t'
)


# Create dictionary that maps ACCESSION.VERSION's to seqIDs
print('Building `acc_seqIDs_dict`')
acc_seqIDs_dict = make_acc_seqIDs_dict(fasta_seqs_fpath, seqkit_fpath)
# accs_with_16S_genes = set(acc_seqIDs_dict.keys())
print('`acc_seqIDs_dict` is built')

# Create tuple of Assembly IDs
ass_ids = tuple(
    set(
        assm_acc_df['ass_id']
    )
)


# == Proceed ==

with open(per_genome_outfpath, 'wt') as per_genome_outfile:

    # Write headers to output files
    per_genome_outfile.write('ass_id\taccs\ttaxID\n')

    # Iterate over Assembly IDs
    for i, ass_id in enumerate(ass_ids):

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        # Get all ACCESSION.VERSION's of current assembly
        accs = tuple(
            assm_acc_df[assm_acc_df['ass_id'] == ass_id]['acc']
        )

        # Try elink 3 times
        error = True
        n_errors = 0
        while error:
            try:
                handle = Entrez.elink(
                    dbfrom='assembly',
                    db='taxonomy',
                    id=ass_id
                )
                records = Entrez.read(handle)
                handle.close()
            except OSError:
                n_errors =+ 1
                if n_errors == 3:
                    records = list()
                    print(f'Oh no, it error...: {ass_id}')
                # end if
            else:
                error = False
            # end try
        # end while

        # Retrieve taxID from response
        try:
            taxID = records[0]['LinkSetDb'][0]['Link'][0]['Id']
        except IndexError as err:
            print(f'Error on {accs}: {err}')
            taxID = 'NA'
        # end try

        # Write to per-genome output file
        per_genome_outfile.write(f'{ass_id}\t{";".join(accs)}\t{taxID}\n')

        # Wait a bit: we don't want NCBI to ban us :)
        time.sleep(0.4)
    # end for
# end with

print('\nCompleted!')
print(per_genome_outfpath)
print(f'\n|=== EXITTING SCRIPT `{__file__}` ===|\n')
