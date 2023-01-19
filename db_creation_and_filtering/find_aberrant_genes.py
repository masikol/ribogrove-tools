#!/usr/bin/env python3

# The script finds aberrant genes: truncated genes and genes with large deletions.

## Command line arguments

### Input files:
# 1. `-f / --fasta-seqs-file` -- an input fasta file of all collected SSU gene sequences.
#   This file is the output of the script `extract_16S.py`. Mandatory.
# 2. `--NNN-fail-seqIDs` -- a file of seqIDs which didn't pass NNN filter, one per line.
#   This file is the output of the script `find_NNN.py`. Mandatory.
# 3. `-t / --cmscan-tblout` -- a TSV file (`.tblout`) of comparison statistics.
#   This file is the output of the script `compare_all_seqs_to_cm.py`. Mandatory.
# 4. `-c / --conserved-regions-fasta` -- a fasta file of 16S rDNA conserved regions
#   from work "How conserved are the conserved 16S-rRNA regions? (https://peerj.com/articles/3036/)"
#   (table 5, "NR" sequences, but flanking Ns should be removed). Optional.

### Output files:
# 1. `-o / --outdir` -- output directory, where all output files will be stored. Mandatory.
# This directory will contain the following files:
# - `pivotal_genes.tsv` -- a TSV file, with seqIDs of pivotal genes and corresponding cmscan scores.
# - `pident_pivotal_genes.tsv` -- a TSV file containing perpents of identity
#   (and some other statistics) of pairwise alignments of pivotal genes vs non-pivotal genes.
# - `insertions.tsv` -- a TSV file containing information about discovered insertions
#   longer than `--deletion-len-threshold`. Quite counterintuitive but anyway.
# - `deletions.tsv` -- a TSV file containing information about discovered deletions
#   longer than `--deletion-len-threshold`.
# - `non_aberrant_seqIDs.txt` -- seqIDs of discovered non-aberrant genes, one per line.
# - `aberrant_seqIDs.txt` -- seqIDs of discovered aberrant genes, one per line.

### Dependencies:
# 1. `--muscle` -- [MUSCLE](https://www.drive5.com/muscle/) aligner executable
#   for multiple sequence alignment. Mandatory.

### Parameters:
# 1. `--deletion-len-threshold` -- a threshold value for length of an deletion.
#   If deletion has length higher than this value, gene having this deletion will
#   be regarded as an aberrant gene. Integer number > 0. Mandatory.

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import re
import sys
import argparse
import operator
from io import StringIO
import subprocess as sp
from functools import reduce
from typing import Sequence, Dict, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord

from read_and_filter_fasta import read_and_filter_fasta


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--fasta-seqs-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '--NNN-fail-seqIDs',
    help='input file of seqIDs which didn\'t pass NNN filter',
    required=True
)

parser.add_argument(
    '-a',
    '--assm-acc-file',
    help='a TSV file of 4 columns: (`ass_id`, `gi_number`, `acc`, `title`).',
    required=True
)

parser.add_argument(
    '-t',
    '--cmscan-tblout',
    help='TSV file (.tblout) outputed by the script `compare_all_seqs_to_cm.py`',
    required=True
)

parser.add_argument(
    '-c',
    '--conserved-regions-fasta',
    help="""fasta file of NR conserved regions from work
    "How conserved are the conserved 16S-rRNA regions?" (table 5, https://peerj.com/articles/3036/)""",
    required=False
)

# Output files

parser.add_argument(
    '-o',
    '--outdir',
    help='output directory',
    required=True
)


# Dependencies

parser.add_argument(
    '--muscle',
    help='MUSCLE aligner executable',
    required=True
)


# Heuristic's params

parser.add_argument(
    '--deletion-len-threshold',
    help="""Threshold value for length of an deletion.
    If deletion has length higher than this value, gene having this deletion
    is regarded as an aberrant gene.""",
    required=True
)


args = parser.parse_args()


# For convenience
fasta_seqs_fpath = os.path.abspath(args.fasta_seqs_file)
NNN_fail_fpath = os.path.abspath(args.NNN_fail_seqIDs)
assm_acc_fpath = os.path.abspath(args.assm_acc_file)
tblout_fpath = os.path.abspath(args.cmscan_tblout)
muscle_fpath = os.path.abspath(args.muscle)
outdpath = os.path.abspath(args.outdir)


try:
    conserved_regions_fpath = os.path.abspath(args.conserved_regions_fasta)
except TypeError:
    conserved_regions_fpath = None
else:
    if not os.path.exists(conserved_regions_fpath):
        print(f'Error: file `{conserved_regions_fpath}` does not exist!')
        sys.exit(1)
    # end if
# end try


# Check (and assign) value of `--deletion-len-threshold`
try:
    deletion_len_threshold = int(args.deletion_len_threshold)
    if deletion_len_threshold < 0:
        raise ValueError
    # end if
except ValueError:
    print('Error: you entered and invalid `--deletion-len-threshold` value.')
    print(f'Your value: `{args.deletion_len_threshold}`')
    print('It must be integer number >= 0')
    sys.exit(1)
# end try


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, assm_acc_fpath, tblout_fpath, muscle_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Check if executables are actually executable
if not os.access(muscle_fpath, os.X_OK):
    print(f'Error: file `{muscle_fpath}` is not executable!')
    sys.exit(1)
# end if

# Create output directory if needed
if not os.path.isdir(outdpath):
    try:
        os.makedirs(outdpath)
    except OSError as err:
        print(f'Error: cannot create directory `{outdpath}`')
        sys.exit(1)
    # end try
# end if

print(fasta_seqs_fpath)
print(NNN_fail_fpath)
print(assm_acc_fpath)
print(tblout_fpath)
print(conserved_regions_fpath)
print(muscle_fpath)
print()


def select_gene_seqs(ass_id: str,
                     seq_records: Sequence[str]) -> Dict[str, SeqRecord]:
    selected_seq_records = filter(
        lambda r: get_ass_id_from_seq_record(r) == ass_id,
        seq_records
    )
    # Make result dictionary and return it
    return {r.id: r for r in selected_seq_records}
# end def


def get_ass_id_from_seq_record(seq_record):
    return int(seq_record.id.partition(':')[0][2:])
# end def


def pairwise_align(
    pivotal_seq_record: SeqRecord,
    seq_record: SeqRecord,
    muscle_fpath: str) -> Tuple[SeqRecord, SeqRecord]:
    # Function performs pairwise global alignment of a pivotal gene (pivotal_seq_record)
    #    and some another gene (seq_record).
    # It does pairwise with MUSCLE. It's quite weird, but we'll do it
    #    for the sake of uniformity.

    # Configure command
    cmd = f'{muscle_fpath} -diags -quiet'
    pipe = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)

    # Configure input fasta string
    fasta_str = f'>{pivotal_seq_record.id}\n{str(pivotal_seq_record.seq)}\n>{seq_record.id}\n{str(seq_record.seq)}\n'

    # Write fasta string to stdin
    pipe.stdin.write(fasta_str.encode('ascii'))

    # Run command
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error aligning gene seqs!')
        print(stdout_stderr[1].decode('utf-8'))
        print(f'pivotal_seq_record.id = {pivotal_seq_record.id}')
        print(f'seq_record.id = {seq_record.id}')
        print(f'\t{cmd}')
        sys.exit(1)
    else:
        # Obtain result alignment
        aln_str = stdout_stderr[0].decode('utf-8')
    # end if

    # Parse result alignment
    aln_io = StringIO(aln_str)
    aln_records = list(SeqIO.parse(aln_io, 'fasta'))
    aln_io.close()

    # Extract aligned strings of pivotal gene and another gene alone
    pivotal_aln_record = next(filter(lambda r: r.id == pivotal_seq_record.id, aln_records))
    aln_record = next(filter(lambda r: r.id != pivotal_seq_record.id, aln_records))

    return pivotal_aln_record, aln_record
# end def


def bases_identical(b1: str, b2: str) -> bool:
    # Function checks if two bases are identical
    return 1 if b1 == b2 else 0
# end def


def pairwise_percent_identity(aln_record_1: SeqRecord, aln_record_2: SeqRecord) -> float:
    # Function calculates percent of identity of a pairwise alignment (pident).

    # Extract aligned sequences as type `str`
    seq_1 = str(aln_record_1.seq)
    seq_2 = str(aln_record_2.seq)

    # Find length of the shortest sequence
    min_len = min(
        map(
            lambda x: len(x),
            [
                seq_1.replace('-', ''),
                seq_2.replace('-', '')
            ]
        )
    )

    # Calculate and return pident
    return sum( (bases_identical(b1, b2) for b1, b2 in zip(seq_1, seq_2)) ) \
           / min_len
# end def


def count_gaps(seq_record: SeqRecord) -> int:
    # Function counts gaps in an anigned sequence
    return str(seq_record.seq).count('-')
# end def



def find_insertions_and_deletions(
    pivotal_aln_record: SeqRecord,
    aln_record: SeqRecord,
    deletion_len_threshold: int) -> Tuple[Tuple[int, int, str], Tuple[int, int]]:
    # Function find long insertions and deletions

    # `deletion_len_threshold` minuses in a row
    indel_pattern = r'[-]{%d,}' % int(deletion_len_threshold+1)

    # = Find insertions =

    insertion_matchobjs = tuple(
        re.finditer(indel_pattern, str(pivotal_aln_record.seq))
    )

    insertions = list()
    for obj in insertion_matchobjs:
        ins_start = obj.start()
        ins_end = obj.end()
        ins_seq = str(aln_record.seq)[ins_start : ins_end] # inserted sequence
        insertions.append(
            # To 1-based, left-closed, right-closed
            (ins_start+1, ins_end, ins_seq)
        )
    # end for

    # = Find deletions =

    deletion_matchobjs = tuple(
        re.finditer(indel_pattern, str(aln_record.seq))
    )

    # To 1-based, left-closed, right-closed
    deletions = [(obj.start()+1, obj.end()) for obj in deletion_matchobjs]

    return insertions, deletions
# end def


def set_acc(row):
    row['acc'] = row['query_name'].split(':')[1]
    return row
# end def


if not conserved_regions_fpath is None:

    def find_conserved_regions(seq_record: SeqRecord,
                               conserved_seq_records: Sequence[SeqRecord]) -> None:
        # Function returns seqIDs of those conserved regions, which are present in
        #   sequence of record `seq_record.

        # Case sequence to type `str`
        seq = str(seq_record.seq)

        # Find conserved regions, which are present in `seq`
        present_conserved_regions = filter(
            # Coordinates of conserved region's occurence(s) are stored in
            #   list returned by SeqUtils.nt_search. If there is no occurence,
            #   this list contain single element -- template sequence (`seq`).
            lambda cons_rec: len(SeqUtils.nt_search(seq, str(cons_rec.seq))) > 1,
            conserved_seq_records
        )

        # Return seqIDs of present conserved regions
        return set(
            map(
                lambda cons_rec: cons_rec.id,
                present_conserved_regions
            )
        )
    # end def
# end if


# Output files
pivotal_genes_fpath = os.path.join(outdpath, 'pivotal_genes.tsv')
pident_outfpath = os.path.join(outdpath, 'pident_pivotal_genes.tsv')
insertions_outfpath = os.path.join(outdpath, 'insertions.tsv')
deletions_outfpath = os.path.join(outdpath, 'deletions.tsv')
aberrant_seqIDs_fpath = os.path.join(outdpath, 'aberrant_seqIDs.txt')



# == Proceed ==

# Read conserved regions' sequences
if not conserved_regions_fpath is None:
    conserved_seq_records = tuple(SeqIO.parse(conserved_regions_fpath, 'fasta'))
# end if

# Read pivotal genes' dataframe
# pivotal_genes_df = pd.read_csv(pivotal_genes_fpath, sep='\t')

# Reade per-replicon genes statistics
assm_acc_df = pd.read_csv(assm_acc_fpath, sep='\t')

# Read cmscan's tblout file
tblout_df = pd.read_csv(tblout_fpath, sep='\t')

tblout_df['acc'] = np.repeat(None, tblout_df.shape[0])
tblout_df = tblout_df.apply(set_acc, axis=1)

tblout_df = tblout_df.merge(
    assm_acc_df[['acc', 'ass_id']],
    on='acc',
    how='left'
)

# Get unique Assembly IDs
ass_ids = tuple(set(assm_acc_df['ass_id']))

# Read input genes sequences
seq_records = read_and_filter_fasta(
    fasta_seqs_fpath,
    filter_fpaths=[NNN_fail_fpath,]
)


with open(pivotal_genes_fpath, 'wt') as pivotal_genes_outfile, \
     open(pident_outfpath, 'wt') as pident_outfile, \
     open(insertions_outfpath, 'wt') as insertions_outfile, \
     open(deletions_outfpath, 'wt') as deletions_outfile, \
     open(aberrant_seqIDs_fpath, 'wt') as aberrant_seqIDs_outfile:

    # Write headers to output files
    pivotal_genes_outfile.write('ass_id\tpivotal_seqID\tmax_score\n')
    pident_outfile.write('ass_id\tpivotal_seqID\tseqID\tpident\tn_insert_bases\tn_delet_bases\tseqID_is_also_pivotal\n')
    insertions_outfile.write('ass_id\tpivotal_seqID\tseqID\tpivotal_gene_len\tgene_len\tstart\tend\tseq\n')
    deletions_outfile.write('ass_id\tpivotal_seqID\tseqID\tpivotal_gene_len\tgene_len\tstart\tend\n')

    # Iterate over assemblies
    for i, ass_id in enumerate(ass_ids):
        print(f'\rDoing genome #{i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        # Select SSU genes from current genome
        selected_seq_records = select_gene_seqs(ass_id, seq_records)

        if len(selected_seq_records) == 0:
            continue
        # end if

        # Get all seqIDs of genes from current genome
        seqIDs = set(selected_seq_records.keys())

        # Select rows corresponding to current assembly
        curr_ass_df = tblout_df[tblout_df['ass_id'] == ass_id]
        # Save seqIDs of aligned sequences into a set for faster search
        aligned_seqIDs = set(curr_ass_df['query_name'])

        # Select a priori aberrant genes: truncated and analigned ones
        apriori_aberrant_seqIDs = set(
            curr_ass_df[curr_ass_df['trunc'] != 'no']['query_name']
        ) | set(
            filter(
                lambda x: not x in aligned_seqIDs,
                seqIDs
            )
        )

        # Remove a priori aberrant genes from `curr_ass_df`
        curr_ass_df = curr_ass_df.query('not query_name in @apriori_aberrant_seqIDs')

        # Select pivotal genes
        max_score = curr_ass_df['score'].max()
        pivotal_seqIDs = set(curr_ass_df[curr_ass_df['score'] > (max_score - 1e-6)]['query_name'])

        # Record seqIDs of pivotal genes
        for pivotal_seqID in pivotal_seqIDs:
            pivotal_genes_outfile.write(f'{ass_id}\t{pivotal_seqID}\t{max_score}\n')
        # end for

        # Here sets of seqIDs of aberrant genes will be stored, in comparison to
        #   each pivotal gene
        aberrant_seqIDs_setlist = list()

        # Iterate over ivotal genes
        for pivotal_seqID in pivotal_seqIDs:
            pivotal_seq_record = selected_seq_records[pivotal_seqID]

            # Here seqIDs of aberrant genes will be stored, but only in comparison to
            #   current pivotal gene
            curr_aberrant_seqIDs = set()

            # Save conserved regions present in current pivotal gene
            if not conserved_regions_fpath is None:
                conserv_IDs_in_pivotal_gene = find_conserved_regions(
                    pivotal_seq_record,
                    conserved_seq_records
                )
            # end if

            # Iterate over all genes i nthe genome, except for current pivotal gene
            for seqID in seqIDs - {pivotal_seqID}:

                seq_record = selected_seq_records[seqID]

                # Perform pairwise alignment of current gene and current pivotal gene
                pivotal_aln_record, aln_record = pairwise_align(pivotal_seq_record, seq_record, muscle_fpath)

                # Save some statistics of performaed alignment
                pident = pairwise_percent_identity(pivotal_aln_record, aln_record)
                n_insert_bases = count_gaps(pivotal_aln_record)
                n_delet_bases = count_gaps(aln_record)
                seqID_is_also_pivotal = seqID in pivotal_seqIDs

                # Write these statistics of pairwise alignment
                pident_outfile.write(f'{ass_id}\t{pivotal_seqID}\t{seqID}\t{pident}\t')
                pident_outfile.write(f'{n_insert_bases}\t{n_delet_bases}\t{1 if seqID_is_also_pivotal else 0}\n')

                # Find long insertions and deletions
                insertions, deletions = find_insertions_and_deletions(
                    pivotal_aln_record,
                    aln_record,
                    deletion_len_threshold
                )

                # Record insertions
                for insertion in insertions:
                    insertions_outfile.write(f'{ass_id}\t{pivotal_seqID}\t{seqID}\t')
                    insertions_outfile.write(f'{len(pivotal_seq_record.seq)}\t{len(seq_record.seq)}\t')
                    insertions_outfile.write(f'{insertion[0]}\t{insertion[1]}\t{insertion[2]}\n')
                # end for

                # Record deletions
                for deletion in deletions:
                    deletions_outfile.write(f'{ass_id}\t{pivotal_seqID}\t{seqID}\t')
                    deletions_outfile.write(f'{len(pivotal_seq_record.seq)}\t{len(seq_record.seq)}\t')
                    deletions_outfile.write(f'{deletion[0]}\t{deletion[1]}\n')
                # end for

                # Save conserved regions present in current gene
                if not conserved_regions_fpath is None:
                    conserv_IDs_in_curr_gene = find_conserved_regions(
                        seq_record,
                        conserved_seq_records
                    )

                    # If this flag is True, some conserved regions are missing, in comparison to
                    #   current pivotal gene
                    missing_conserved_regions = len(conserv_IDs_in_curr_gene) < len(conserv_IDs_in_pivotal_gene)
                else:
                    missing_conserved_regions = False
                # end if

                # If there are long indels or some conserved regions are missing, current gene
                #   is aberrant in comparison to current pivotal gene
                # if len(insertions) != 0 or len(deletions) != 0 or missing_conserved_regions:
                # if len(deletions) != 0:
                # if missing_conserved_regions:
                if len(deletions) != 0 or missing_conserved_regions:
                    curr_aberrant_seqIDs.add(seqID)
                # end if
            # end for

            # Add set of current aberrant genes to our list
            aberrant_seqIDs_setlist.append(curr_aberrant_seqIDs)
        # end for

        # Aberrant genes will be those genes, which are aberrant in comparison to
        #   all pivotal genes
        final_aberrant_seqIDs = apriori_aberrant_seqIDs
        if len(aberrant_seqIDs_setlist) != 0:
            final_aberrant_seqIDs = final_aberrant_seqIDs | reduce(operator.and_, aberrant_seqIDs_setlist)
        # end if

        non_aberrant_seqIDs = seqIDs - final_aberrant_seqIDs

        # Record seqIDs of aberrant genes
        for seqID in final_aberrant_seqIDs:
            aberrant_seqIDs_outfile.write(f'{seqID}\n')
        # end for
    # end for
# end with

print('\nCompleted!')
print(pivotal_genes_fpath)
print(pident_outfpath)
print(insertions_outfpath)
print(deletions_outfpath)
print(aberrant_seqIDs_fpath)

print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
