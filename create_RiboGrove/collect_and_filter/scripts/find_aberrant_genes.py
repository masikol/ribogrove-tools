#!/usr/bin/env python3

# The script finds aberrant genes: truncated genes and genes with large deletions.

## Command line arguments

### Input files:
# 1. `-f / --fasta-seqs-file` -- an input fasta file of all collected SSU gene sequences.
#   This file is the output of the script `extract_16S.py`. Mandatory.
# 2. `--ribotyper-fail-seqIDs` -- a file outputted by the script `find_ribotyper_fail_seqs.py`.
#   It lists seqIDs of sequences which failed the ribotyper filter, one per line.
#   Mandatory.
# 3. `-a / --in-asm-sum` -- an assembly summary file after the 2nd step of filtering
#   Mandatory.
# 4. `-l / --ribotyper-long-out-tsv` -- input .long.out.tsv file
#   outputted by the script `check_seqs_with_ribotyper.py`.
#   Mandatory.

### Output files:
# 1. `-o / --outdir` -- output directory, where all output files will be stored.
# This directory will contain the following files:
# - `pivotal_genes.tsv` -- a TSV file, with seqIDs of pivotal genes and corresponding cmscan scores.
# - `pident_pivotal_genes.tsv` -- a TSV file containing perpents of identity
#   (and some other statistics) of pairwise alignments of pivotal genes vs non-pivotal genes.
# - `insertions.tsv` -- a TSV file containing information about discovered insertions
#   longer than `--deletion-len-threshold`. Quite counterintuitive but anyway.
# - `deletions.tsv` -- a TSV file containing information about discovered deletions
#   longer than `--deletion-len-threshold`.
# - `aberrant_seqIDs.txt` -- seqIDs of discovered aberrant genes, one per line.
#   Mandatory.

### Dependencies:
# 1. `--muscle` -- [MUSCLE](https://www.drive5.com/muscle/) aligner executable
#   for multiple sequence alignment. Mandatory.

### Parameters:
# 1. `--deletion-len-threshold` -- a threshold value for length of an deletion.
#   If deletion has length higher than this value, gene having this deletion will
#   be regarded as an aberrant gene. Integer number > 0. Mandatory.

### "Cached" files:
# 1. `--prev-final-fasta` -- a final fasta file of 16S rRNA gene sequences
#   from the previous RiboGrove release.
#   Optional.
# 2. `--prev-aberrant-seqIDs` -- file `aberrant_seqIDs.txt`
#   from the previous RiboGrove release.
#   Optional.


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
    '--fasta-seqs-file',
    help='fasta file of 16S rRNA gene sequences',
    required=True
)

parser.add_argument(
    '--ribotyper-fail-seqIDs',
    help="""a file outputted by the script `find_ribotyper_fail_seqs.py`.
   It lists seqIDs of sequences which failed the ribotyper filter, one per line.""",
    required=True
)

parser.add_argument(
    '-a',
    '--in-asm-sum',
    help='an assembly summary file after the 2nd step of filtering',
    required=True
)

parser.add_argument(
    '-l',
    '--ribotyper-long-out-tsv',
    help="""input .long.out.tsv file
    outputted by the script `check_seqs_with_ribotyper.py`""",
    required=True
)

# Cache

parser.add_argument(
    '--prev-final-fasta',
    help="""a final fasta file of 16S rRNA gene sequences
    from the previous RiboGrove release""",
    required=False
)

parser.add_argument(
    '--prev-aberrant-seqIDs',
    help="""file `aberrant_seqIDs.txt`
    from the previous RiboGrove release""",
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


# == Import them now ==
import re
import sys
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

import src.rg_tools_IO as rgIO
from src.ribogrove_seqID import parse_asm_acc


# For convenience
fasta_seqs_fpath = os.path.abspath(args.fasta_seqs_file)
ribotyper_fail_fpath = os.path.abspath(args.ribotyper_fail_seqIDs)
asm_sum_fpath = os.path.abspath(args.in_asm_sum)
ribotyper_long_fpath = os.path.abspath(args.ribotyper_long_out_tsv)
muscle_fpath = os.path.abspath(args.muscle)
if not args.prev_final_fasta is None \
   and not args.prev_aberrant_seqIDs is None:
    prev_final_fasta_fpath = os.path.abspath(args.prev_final_fasta)
    prev_aberrant_seqIDs_fpath = os.path.abspath(args.prev_aberrant_seqIDs)
else:
    prev_final_fasta_fpath = None
    prev_aberrant_seqIDs_fpath = None
# end if
outdpath = os.path.abspath(args.outdir)


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
for fpath in (fasta_seqs_fpath, asm_sum_fpath, ribotyper_long_fpath, muscle_fpath):
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

# Check if cache files are specified
if not prev_final_fasta_fpath is None \
   and not prev_aberrant_seqIDs_fpath is None:
    # A flag variable for convenience
    cache_mode = True
    for f in (prev_final_fasta_fpath, prev_aberrant_seqIDs_fpath):
        if not os.path.exists(f):
            print(f'Error: file `{f}` does not exist')
            sys.exit(1)
        # end if
    # end for
else:
    cache_mode = False
# end if

print(fasta_seqs_fpath)
print(ribotyper_fail_fpath)
print(asm_sum_fpath)
print(ribotyper_long_fpath)
print(muscle_fpath)
print('deletion_len_threshold = {}'.format(deletion_len_threshold))
if cache_mode:
    print(prev_final_fasta_fpath)
    print(prev_aberrant_seqIDs_fpath)
# end if
print()


def get_prev_seqIDs(prev_aberrant_seqIDs_fpath):
    with open(prev_aberrant_seqIDs_fpath, 'rt') as seqIDs_file:
        lines = seqIDs_file.readlines()
    # end with

    return set(
        map(str.strip, lines)
    )
# end def

def get_cached_aberrant_seqIDs(seq_records, prev_aberrant_seqIDs):
    all_curr_seqIDs = set(
        map(lambda r: r.id, seq_records)
    )
    return all_curr_seqIDs & prev_aberrant_seqIDs
# end def

def make_cached_asm_accs(cached_aberrant_seqIDs, prev_final_fasta_fpath, all_asm_accs):
    cached_aberrant_asm_accs = set(
        map(parse_asm_acc, cached_aberrant_seqIDs)
    )
    prev_final_asm_accs = set(
        map(
            lambda record: parse_asm_acc(record.id),
            SeqIO.parse(prev_final_fasta_fpath, 'fasta')
        )
    )
    return all_asm_accs & (cached_aberrant_asm_accs | prev_final_asm_accs)
# end def




def select_gene_seqs(asm_acc: str,
                     seq_records: Sequence[str]) -> Dict[str, SeqRecord]:
    selected_seq_records = filter(
        lambda r: get_asm_acc_from_seq_record(r) == asm_acc,
        seq_records
    )
    # Make result dictionary and return it
    return {r.id: r for r in selected_seq_records}
# end def


def get_asm_acc_from_seq_record(seq_record):
    return parse_asm_acc(seq_record.id)
# end def


def pairwise_align(pivotal_seq_record: SeqRecord,
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


def pairwise_percent_identity(aln_record_1: SeqRecord,
                              aln_record_2: SeqRecord) -> float:
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



def find_insertions_and_deletions(pivotal_aln_record: SeqRecord,
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


def read_rt_long_df(ribotyper_long_fpath, ribotyper_fail_fpath):
    rt_long_df = pd.read_csv(
        ribotyper_long_fpath,
        sep='\t',
        dtype={
            'target': str,
            'pass_fail': str,
            'length': pd.Int32Dtype(),
            'fm': str,
            'fam': str,
            'domain': str,
            'model': str,
            'strnd': str,
            'ht': str,
            'tscore': str,
            'bscore': str,
            's_per_nt': str,
            'bevalue': str,
            'tcov': str,
            'bcov': str,
            'bfrom': str,
            'bto': str,
            'mfrom': str,
            'mto': str,
            'scdiff': str,
            'scd_per_nt': str,
            'model': str,
            'tscore': str,
            'unexpected_features': str,
        }
    )

    # Remove sequences which failed ribotyper filter
    ribotyper_fail_seqIDs = read_ribotyper_failed_seqIDs(ribotyper_fail_fpath)
    rt_long_df = rt_long_df.query('not target in @ribotyper_fail_seqIDs').copy()

    # for NoHist, whose `bscore` is '-'
    rt_long_df['bscore'] = rt_long_df['bscore'].replace({'-': 0.0})
    rt_long_df['bscore'] = rt_long_df['bscore'].map(float)

    rt_long_df['asm_acc'] = np.repeat('', rt_long_df.shape[0])
    rt_long_df = rt_long_df.apply(set_asm_acc, axis=1)
    return rt_long_df
# end def

def read_ribotyper_failed_seqIDs(ribotyper_fail_fpath):
    with open(ribotyper_fail_fpath, 'rt') as infile:
        fail_seqIDs = set(
            map(
                str.strip,
                infile.readlines()
            )
        )
    # end with
    return fail_seqIDs
# end def

def set_asm_acc(row):
    row['asm_acc'] = parse_asm_acc(row['target'])
    return row
# end def


# Output files
pivotal_genes_fpath   = os.path.join(outdpath, 'pivotal_genes.tsv')
pident_outfpath       = os.path.join(outdpath, 'pident_pivotal_genes.tsv')
insertions_outfpath   = os.path.join(outdpath, 'insertions.tsv')
deletions_outfpath    = os.path.join(outdpath, 'deletions.tsv')
aberrant_seqIDs_fpath = os.path.join(outdpath, 'aberrant_seqIDs.txt')



# == Proceed ==

# Reade per-replicon genes statistics
asm_sum_df = rgIO.read_ass_sum_file(asm_sum_fpath)

# Read ribotyper's long.out.tsv file
rt_long_df = read_rt_long_df(ribotyper_long_fpath, ribotyper_fail_fpath)


# Get unique Assembly accessions
asm_accs = set(asm_sum_df['asm_acc'])

# Read input genes sequences
seq_records = rgIO.read_and_filter_fasta(
    fasta_seqs_fpath,
    filter_fpaths=[ribotyper_fail_fpath,]
)

if cache_mode:
    prev_aberrant_seqIDs = get_prev_seqIDs(prev_aberrant_seqIDs_fpath)
    cached_aberrant_seqIDs = get_cached_aberrant_seqIDs(seq_records, prev_aberrant_seqIDs)
    cached_asm_accs = make_cached_asm_accs(cached_aberrant_seqIDs, prev_final_fasta_fpath, asm_accs)
    print('{}/{} genomes are cached'.format(len(cached_asm_accs), len(asm_accs)))
    print('{} genomes left to_process'.format(len(asm_accs - cached_asm_accs)))
    del prev_aberrant_seqIDs
else:
    cached_aberrant_seqIDs = set()
    cached_asm_accs = set()
# end if


with open(pivotal_genes_fpath, 'wt') as pivotal_genes_outfile, \
     open(pident_outfpath, 'wt') as pident_outfile, \
     open(insertions_outfpath, 'wt') as insertions_outfile, \
     open(deletions_outfpath, 'wt') as deletions_outfile, \
     open(aberrant_seqIDs_fpath, 'wt') as aberrant_seqIDs_outfile:

    # Write headers to output files
    pivotal_genes_outfile.write('asm_acc\tpivotal_seqID\tmax_score\n')
    pident_outfile.write('asm_acc\tpivotal_seqID\tseqID\tpident\tn_insert_bases\tn_delet_bases\tseqID_is_also_pivotal\n')
    insertions_outfile.write('asm_acc\tpivotal_seqID\tseqID\tpivotal_gene_len\tgene_len\tstart\tend\tseq\n')
    deletions_outfile.write('asm_acc\tpivotal_seqID\tseqID\tpivotal_gene_len\tgene_len\tstart\tend\n')

    # Iterate over assemblies
    asm_accs_left = asm_accs - cached_asm_accs
    for i, asm_acc in enumerate(asm_accs_left):
        print(
            '\r{} -- Doing genome #{}/{}: {}' \
                .format(get_time(), i+1, len(asm_accs_left), asm_acc),
            end=' '*10
        )

        # Select SSU genes from current genome
        selected_seq_records = select_gene_seqs(asm_acc, seq_records)

        if len(selected_seq_records) == 0:
            continue
        # end if

        # Get all seqIDs of genes from current genome
        curr_seqIDs = set(selected_seq_records.keys())

        # Select rows corresponding to current assembly
        curr_ass_df = rt_long_df[rt_long_df['asm_acc'] == asm_acc]

        # Select pivotal genes
        max_score = curr_ass_df['bscore'].max()
        pivotal_seqIDs = set(
            curr_ass_df[curr_ass_df['bscore'] > (max_score - 1e-6)]['target']
        )

        # Record seqIDs of pivotal genes
        for pivotal_seqID in pivotal_seqIDs:
            pivotal_genes_outfile.write(f'{asm_acc}\t{pivotal_seqID}\t{max_score}\n')
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

            # Iterate over all genes i nthe genome, except for current pivotal gene
            for seqID in curr_seqIDs - {pivotal_seqID}:

                seq_record = selected_seq_records[seqID]

                # Perform pairwise alignment of current gene and current pivotal gene
                pivotal_aln_record, aln_record = pairwise_align(pivotal_seq_record, seq_record, muscle_fpath)

                # Save some statistics of performaed alignment
                pident = pairwise_percent_identity(pivotal_aln_record, aln_record)
                n_insert_bases = count_gaps(pivotal_aln_record)
                n_delet_bases = count_gaps(aln_record)
                seqID_is_also_pivotal = seqID in pivotal_seqIDs

                # Write these statistics of pairwise alignment
                pident_outfile.write(f'{asm_acc}\t{pivotal_seqID}\t{seqID}\t{pident}\t')
                pident_outfile.write(f'{n_insert_bases}\t{n_delet_bases}\t{1 if seqID_is_also_pivotal else 0}\n')

                # Find long insertions and deletions
                insertions, deletions = find_insertions_and_deletions(
                    pivotal_aln_record,
                    aln_record,
                    deletion_len_threshold
                )

                # Record insertions
                for insertion in insertions:
                    insertions_outfile.write(f'{asm_acc}\t{pivotal_seqID}\t{seqID}\t')
                    insertions_outfile.write(f'{len(pivotal_seq_record.seq)}\t{len(seq_record.seq)}\t')
                    insertions_outfile.write(f'{insertion[0]}\t{insertion[1]}\t{insertion[2]}\n')
                # end for

                # Record deletions
                for deletion in deletions:
                    deletions_outfile.write(f'{asm_acc}\t{pivotal_seqID}\t{seqID}\t')
                    deletions_outfile.write(f'{len(pivotal_seq_record.seq)}\t{len(seq_record.seq)}\t')
                    deletions_outfile.write(f'{deletion[0]}\t{deletion[1]}\n')
                # end for

                # If there are long indels, current gene
                #   is aberrant in comparison to current pivotal gene
                if len(deletions) != 0:
                    curr_aberrant_seqIDs.add(seqID)
                # end if
            # end for

            # Add set of current aberrant genes to our list
            aberrant_seqIDs_setlist.append(curr_aberrant_seqIDs)
        # end for

        # Aberrant genes will be those genes, which are aberrant in comparison to
        #   all pivotal genes
        final_aberrant_seqIDs = set()
        if len(aberrant_seqIDs_setlist) != 0:
            final_aberrant_seqIDs = final_aberrant_seqIDs | reduce(operator.and_, aberrant_seqIDs_setlist)
        # end if

        # Record seqIDs of aberrant genes
        for seqID in final_aberrant_seqIDs:
            aberrant_seqIDs_outfile.write(f'{seqID}\n')
        # end for
    # end for

    if cache_mode:
        aberrant_seqIDs_outfile.write('\n'.join(cached_aberrant_seqIDs))
        aberrant_seqIDs_outfile.write('\n')
    # end if
# end with

print('\nCompleted!')
print(pivotal_genes_fpath)
print(pident_outfpath)
print(insertions_outfpath)
print(deletions_outfpath)
print(aberrant_seqIDs_fpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
