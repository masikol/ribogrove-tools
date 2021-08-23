#!/usr/bin/env python3

# Script does multiple things:
# 1. It finds aberrant genes: truncated genes and genes, which diffre from pivotal genes greatly.
# 2. It records all heterogeneity (indels, per-base intragenomic entropy, percents of identity)

# Input files:
# 1. Fasta file of genes sequences (-f/--fasta-seqs-file).
# 2. TSV file if per-replicon genes statistics (-s/--genes-stats-file).
# 3. TSV file TSV file (with header) info about pivotal genes
#    (it is output file of script `find_pivotal_gens.py`).
# 4. Fasta file of NR conserved regions from work
#    "How conserved are the conserved 16S-rRNA regions?"
#    (table 5, https://peerj.com/articles/3036/)
#    -c/--conserved-regions-fasta

# Output files (they all will be stored in output directory):
# 1. pident_pivotal_genes.tsv -- TSV file containing pidents (and some other statistics)
#    of pairwise alignments of pivotal genes vs non-pivotal genes.
# 2. insertions.tsv -- TSV file containing information about discovered insertions
#    longer than `--indel-len-threshold`.
# 3. deletions.tsv -- TSV file containing information about discovered deletions
#    longer than `--indel-len-threshold`.
# 4. entropy.tsv -- TSV file containing per-position intragenomic entropy,
#    excluding aberrant genes.
# 5. aberrant_seqIDs.txt -- seqIDs of aberrant genes, one per line.

# Dependencies:
# 1. MUSCLE aligner (--muscle).

# Parameters:
# 1. Threshold value for length of an indel (--indel-len-threshold).
#    If indel has length higher than this value, gene having this indel
#    is regarded as an aberrant gene.


import os
import re
import sys
import math
import argparse
import operator
from io import StringIO
import subprocess as sp
from array import array
from operator import add
from functools import reduce
from typing import Sequence, Dict, Tuple, List, TextIO

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils
from Bio.SeqRecord import SeqRecord



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
    '-s',
    '--genes-stats-file',
    help='TSV file (with header) containing per-replicons SSU gene statistics',
    required=True
)

parser.add_argument(
    '-p',
    '--pivotal-genes-file',
    help="""TSV file (with header) info about pivotal genes
    (it is output file of script `find_pivotal_gens.py`)""",
    required=True
)

parser.add_argument(
    '-c',
    '--conserved-regions-fasta',
    help="""fasta file of NR conserved regions from work
    "How conserved are the conserved 16S-rRNA regions?" (table 5, https://peerj.com/articles/3036/)""",
    required=True
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
    help='muscle executable',
    required=True
)


# Heuristic's params

parser.add_argument(
    '--indel-len-threshold',
    help="""  """,
    required=True
)


args = parser.parse_args()


# For convenience
fasta_seqs_fpath = os.path.abspath(args.fasta_seqs_file)
genes_stats_fpath = os.path.abspath(args.genes_stats_file)
pivotal_genes_fpath = os.path.abspath(args.pivotal_genes_file)
conserved_regions_fpath = os.path.abspath(args.conserved_regions_fasta)
muscle_fpath = os.path.abspath(args.muscle)
outdpath = os.path.abspath(args.outdir)


# Check (and assign) value of `--indel-len-threshold`
try:
    indel_len_threshold = int(args.indel_len_threshold)
    if indel_len_threshold < 0:
        raise ValueError
    # end if
except ValueError:
    print('Error: you entered and invalid `--indel-len-threshold` value.')
    print(f'Your value: `{args.indel_len_threshold}`')
    print('It must be integer number >= 0')
    sys.exit(1)
# end try


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, genes_stats_fpath, pivotal_genes_fpath, conserved_regions_fpath, muscle_fpath):
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
print(genes_stats_fpath)
print(pivotal_genes_fpath)
print(conserved_regions_fpath)
print(muscle_fpath)
print(outdpath)
print(indel_len_threshold)

# Output files
pident_outfpath = os.path.join(outdpath, 'pident_pivotal_genes.tsv')
insertions_outfpath = os.path.join(outdpath, 'insertions.tsv')
deletions_outfpath = os.path.join(outdpath, 'deletions.tsv')
entropy_outfpath = os.path.join(outdpath, 'entropy.tsv')
aberrant_seqIDs_fpath = os.path.join(outdpath, 'aberrant_seqIDs.txt')



sys.exit(1)


pivotal_genes_fpath = '/mnt/1.5_drive_0/16S_scrubbling/pivotal_genes.tsv'
genes_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_stats_no_NN.tsv'
fasta_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'
conserved_regions_fpath = '/mnt/1.5_drive_0/16S_scrubbling/consensus_seqs/conserved_regions_NR.fasta'

muscle_fpath = '/home/cager/Misc_soft/muscle3.8.31'

pident_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/pident_pivotal_genes.tsv'
insertions_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/insertions.tsv'
deletions_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/deletions.tsv'
entropy_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/entropy.tsv'
aberrant_seqIDs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/aberrant_seqIDs.txt'


def select_gene_seqs(ass_id: str,
    seq_records: Sequence[str],
    stats_df: pd.DataFrame) -> Dict[str, SeqRecord]:

    # Get ACCESSION.VERSION's for current assembly
    accs = set(stats_df[stats_df['ass_id'] == ass_id]['acc'])

    # Filter genes from current genome
    selected_seq_records = tuple(
        filter(
            lambda r: r.id.partition(':')[0] in accs,
            seq_records
        )
    )

    # Make result dictionary and return it
    return {r.id: r for r in selected_seq_records}
# end def select_gene_seqs


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
# end def pairwise_align


def bases_identical(b1: str, b2: str) -> bool:
    # Function checks if two bases are identical
    return 1 if b1 == b2 else 0
# end def bases_identical


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
# end def pairwise_percent_identity


def count_gaps(seq_record: SeqRecord) -> int:
    # Function counts gaps in an anigned sequence
    return str(seq_record.seq).count('-')
# end def count_gaps


def do_msa(seq_records: Sequence[SeqRecord], muscle_fpath: str) -> List[SeqRecord]:
    # Function does Multiple Sequence Alignment

    # Configure command
    cmd = f'{muscle_fpath} -quiet -diags'

    # Configure input fasta string for MSA
    fasta_str = '\n'.join(
        (f'>{r.id}\n{str(r.seq)}' for r in seq_records)
    ) + '\n'

    # Create pipe
    pipe = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    # Write inpit fasta record to stdin
    pipe.stdin.write(fasta_str.encode('ascii'))

    # Run command
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error while doing msa!')
        print('seqIDs:')
        print(' '.join( [r.id for r in seq_records] ))
        print(stdout_stderr[1].decode('utf-8'))
    else:
        # Get StringIO of result alignment
        msa_io = StringIO(stdout_stderr[0].decode('utf-8'))
    # end if

    # Parse alignment
    msa_records = list(SeqIO.parse(msa_io, 'fasta'))
    msa_io.close()

    return msa_records
# end def do_msa


def get_aln_column(i: int, seqs: Sequence[str]) -> Sequence[str]:
    # Function extracts single column from multiple sequence alignment
    return tuple(map(lambda s: s[i], seqs))
# end


def calc_entropy(msa_records: Sequence[SeqRecord]) -> array[float]:
    # Function calculates intragenomic per-base entropy.

    n_seqs = len(msa_records)

    # Extract aligned sequences as type `str`
    seqs = tuple(map(lambda x: str(x.seq), msa_records))
    seq_length = len(seqs[0]) # they are aligned and have the same length

    # Init entropy array
    entropy_arr = array('d', np.repeat(np.nan, seq_length))
    aln_pos = 0

    # Iterate over positions in alignment
    for i in range(seq_length):

        # Extract alignment column of index i
        aln_column = get_aln_column(i, seqs)

        # Count frequences aof all bases at this position
        freqs_arr = tuple(
            map(
                lambda base: aln_column.count(base) / n_seqs,
                set(aln_column)
            )
        )

        # Calculate entropy
        # abs instead of minus in order not to allow "-0.0" values
        entropy_arr[aln_pos] = abs(
            reduce(
                add,
                (freq * math.log(freq, 2) for freq in freqs_arr)
            )
        )
        aln_pos += 1
    # end for

    return entropy_arr
# end def calc_entropy


def calc_and_write_entropy(
    seq_records_for_msa: Sequence[SeqRecord],
    muscle_fpath: str,
    entropy_outfile: TextIO,
    ass_id) -> None:
    # Function-wrapper for calculating intragenomic per-position entropy:
    #   it performs MSA and calculates entropy.

    if len(seq_records_for_msa) > 1:

        # Perform MSA
        msa_records = do_msa(seq_records_for_msa, muscle_fpath)

        # Calculate entropy
        entropy_arr = calc_entropy(msa_records)

        # Write per-base entropy
        for i, entropy in enumerate(entropy_arr):
            entropy_outfile.write(f'{ass_id}\t{i}\t{entropy}\n')
        # end for
    # end if
# end def get_entropy_array


def find_insertions_and_deletions(
    pivotal_aln_record: SeqRecord,
    aln_record: SeqRecord,
    indel_len_threshold: int) -> Tuple[Tuple[int, int, seq], Tuple[int, int]]:
    # Function find long insertions and deletions

    indel_pattern = r'[-]{%d,}' % int(indel_len_threshold+1)

    # = Find insertions =

    insertion_matchobjs = tuple(
        re.finditer(indel_pattern, str(pivotal_aln_record.seq))
    )

    insertions = list()
    for obj in insertion_matchobjs:
        ins_start = obj.start()
        ins_end = obj.end()
        ins_seq = str(aln_record.seq)[ins_start : ins_end]
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
# end def find_insertions_and_deletions


def find_conserved_regions(seq_record, conserved_seq_records):

    seq = str(seq_record.seq)

    return set(
        map(
            lambda cons_rec: cons_rec.id,
            filter(
                lambda cons_rec: len(SeqUtils.nt_search(seq, str(cons_rec.seq))) > 1,
                conserved_seq_records
            )
        )
    )
# end def find_conserved_regions


conserved_seq_records = tuple(SeqIO.parse(conserved_regions_fpath, 'fasta'))


pivotal_genes_df = pd.read_csv(pivotal_genes_fpath, sep='\t')
stats_df = pd.read_csv(genes_stats_fpath, sep='\t')


ass_ids = tuple(set(stats_df['ass_id']))
# ass_ids = [
#     1691841,
    # 131461,
    # 9961891,
    # 3810951,
    # 1442141,
    # 5131611,
    # 9310321,
    # 7359321,
    # 1005941,
    # 5394991,
    # 1491951,
    # 9924121,
# ]

seq_records = tuple(SeqIO.parse(fasta_seqs_fpath, 'fasta'))


with open(pident_outfpath, 'wt') as pident_outfile, \
     open(insertions_outfpath, 'wt') as insertions_outfile, \
     open(deletions_outfpath, 'wt') as deletions_outfile, \
     open(entropy_outfpath, 'wt') as entropy_outfile, \
     open(aberrant_seqIDs_fpath, 'wt') as aberrant_seqIDs_outfile:

    pident_outfile.write('ass_id\tpivotal_seqID\tseqID\tpident\tn_insert_bases\tn_delet_bases\tseqID_is_also_pivotal\n')
    insertions_outfile.write('ass_id\tpivotal_seqID\tseqID\tpivotal_gene_len\tgene_len\tstart\tend\tseq\n')
    deletions_outfile.write('ass_id\tpivotal_seqID\tseqID\tpivotal_gene_len\tgene_len\tstart\tend\n')
    entropy_outfile.write('ass_id\tpivotal_seqID\tpos\tentropy\n')

    for i, ass_id in enumerate(ass_ids):
        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        curr_pivotal_df = pivotal_genes_df[pivotal_genes_df['ass_id'] == ass_id]
        pivotal_gene_num = curr_pivotal_df[~ pd.isnull(curr_pivotal_df['pivotal_gene_seqID'])].shape[0]

        selected_seq_records = select_gene_seqs(ass_id, seq_records, stats_df)

        if pivotal_gene_num != 0:

            seqIDs = set(selected_seq_records.keys())
            pivotal_seqIDs = tuple(curr_pivotal_df['pivotal_gene_seqID'])

            aberrant_seqIDs_setlist = list()

            for pivotal_seqID in pivotal_seqIDs:

                pivotal_seq_record = selected_seq_records[pivotal_seqID]
                curr_aberrant_seqIDs = set()

                conserv_IDs_in_pivotal_gene = find_conserved_regions(
                    pivotal_seq_record,
                    conserved_seq_records
                )

                # print(f'\nPivotal:')
                # print(pivotal_seq_record.id)
                # print(conserv_IDs_in_pivotal_gene)
                # input()

                for seqID in seqIDs - {pivotal_seqID}:

                    seq_record = selected_seq_records[seqID]
                    pivotal_aln_record, aln_record = pairwise_align(pivotal_seq_record, seq_record, muscle_fpath)

                    pident = pairwise_percent_identity(pivotal_aln_record, aln_record)
                    n_insert_bases = count_gaps(pivotal_aln_record)
                    n_delet_bases = count_gaps(aln_record)

                    seqID_is_also_pivotal = seqID in pivotal_seqIDs
                    pident_outfile.write(f'{ass_id}\t{pivotal_seqID}\t{seqID}\t{pident}\t')
                    pident_outfile.write(f'{n_insert_bases}\t{n_delet_bases}\t{1 if seqID_is_also_pivotal else 0}\n')

                    insertions, deletions = find_insertions_and_deletions(pivotal_aln_record, aln_record)

                    for insertion in insertions:
                        insertions_outfile.write(f'{ass_id}\t{pivotal_seqID}\t{seqID}\t')
                        insertions_outfile.write(f'{len(pivotal_seq_record.seq)}\t{len(seq_record.seq)}\t')
                        insertions_outfile.write(f'{insertion[0]}\t{insertion[1]}\t{insertion[2]}\n')
                    # end for

                    for deletion in deletions:
                        deletions_outfile.write(f'{ass_id}\t{pivotal_seqID}\t{seqID}\t')
                        deletions_outfile.write(f'{len(pivotal_seq_record.seq)}\t{len(seq_record.seq)}\t')
                        deletions_outfile.write(f'{deletion[0]}\t{deletion[1]}\n')
                    # end for

                    conserv_IDs_in_curr_gene = find_conserved_regions(
                        seq_record,
                        conserved_seq_records
                    )

                    # print(seq_record.id)
                    # print(conserv_IDs_in_curr_gene)
                    # input()

                    missing_conserved_regions = len(conserv_IDs_in_curr_gene) < len(conserv_IDs_in_pivotal_gene)

                    if len(insertions) != 0 or len(deletions) != 0 or missing_conserved_regions:
                        curr_aberrant_seqIDs.add(seqID)
                    # end if
                # end for

                aberrant_seqIDs_setlist.append(curr_aberrant_seqIDs)
            # end for

            aberrant_seqIDs = reduce(operator.and_, aberrant_seqIDs_setlist)
            for seqID in aberrant_seqIDs:
                aberrant_seqIDs_outfile.write(f'{seqID}\n')
            # end for

            seqIDs_for_msa = seqIDs - aberrant_seqIDs
            seq_records_for_msa = [selected_seq_records[seqID] for seqID in seqIDs_for_msa]
            calc_and_write_entropy(seq_records_for_msa, muscle_fpath, entropy_outfile, ass_id)
        else:
            if tuple(curr_pivotal_df['all_truncated'])[0] == 1:
                aberrant_seqIDs = set(selected_seq_records.keys())
                for seqID in aberrant_seqIDs:
                    aberrant_seqIDs_outfile.write(f'{seqID}\n')
                # end for
            else:
                seq_records_for_msa = selected_seq_records.values()
                calc_and_write_entropy(seq_records_for_msa, muscle_fpath, entropy_outfile, ass_id)
            # end if
        # end if
    # end for
# end with

print('\nCompleted!')
print(pident_outfpath)
print(insertions_outfpath)
print(deletions_outfpath)
print(entropy_outfpath)
print(aberrant_seqIDs_fpath)
