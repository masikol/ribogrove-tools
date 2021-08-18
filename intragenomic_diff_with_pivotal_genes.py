#!/usr/bin/env python3

import os
import re
import sys
import math
import operator
from io import StringIO
import subprocess as sp
from array import array
from operator import add
from functools import reduce

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio import SeqUtils


pivotal_genes_fpath = '/mnt/1.5_drive_0/16S_scrubbling/pivotal_genes.tsv'
stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_stats_no_NN.tsv'
gene_seqs_fasta_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/gene_seqs_no_NN.fasta'
conserved_regions_fpath = '/mnt/1.5_drive_0/16S_scrubbling/consensus_seqs/conserved_regions_NR.fasta'

muscle = '/home/cager/Misc_soft/muscle3.8.31'

pident_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/pident_pivotal_genes.tsv'
insertions_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/insertions.tsv'
deletions_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/deletions.tsv'
entropy_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/entropy.tsv'
aberrant_seqIDs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/aberrant_seqIDs.txt'

# pident_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/test_pident_pivotal_genes.tsv'
# insertions_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/test_insertions.tsv'
# deletions_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/test_deletions.tsv'
# entropy_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/test_entropy.tsv'
# aberrant_seqIDs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/aberrations_and_heterogeneity/test_aberrant_seqIDs.txt'


def select_gene_seqs(ass_id, seq_records, stats_df):

    accs = set(stats_df[stats_df['ass_id'] == ass_id]['acc'])

    selected_seq_records = tuple(
        filter(
            lambda r: r.id.partition(':')[0] in accs,
            seq_records
        )
    )

    return {r.id: r for r in selected_seq_records}
# end def select_gene_seqs


def pairwise_align(pivotal_seq_record, seq_record, muscle):

    cmd = f'{muscle} -diags -quiet'
    pipe = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)

    fasta_str = f'>{pivotal_seq_record.id}\n{str(pivotal_seq_record.seq)}\n>{seq_record.id}\n{str(seq_record.seq)}\n'
    pipe.stdin.write(fasta_str.encode('ascii'))
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error aligning gene seqs!')
        print(stdout_stderr[1].decode('utf-8'))
        print(f'pivotal_seq_record.id = {pivotal_seq_record.id}')
        print(f'seq_record.id = {seq_record.id}')
        print(f'\t{cmd}')
        sys.exit(1)
    else:
        aln_str = stdout_stderr[0].decode('utf-8')
    # end if

    # print(aln_str)

    aln_io = StringIO(aln_str)
    aln_records = list(SeqIO.parse(aln_io, 'fasta'))
    aln_io.close()

    # print(aln_records)

    pivotal_aln_record = next(filter(lambda r: r.id == pivotal_seq_record.id, aln_records))
    aln_record = next(filter(lambda r: r.id != pivotal_seq_record.id, aln_records))

    return pivotal_aln_record, aln_record
# end def pairwise_align


def bases_identical(b1, b2):
    return 1 if b1 == b2 else 0
# end def bases_identical


def pairwise_percent_identity(aln_record_1, aln_record_2):
    seq_1 = str(aln_record_1.seq)
    seq_2 = str(aln_record_2.seq)

    min_len = min(
        map(
            lambda x: len(x),
            [
                seq_1.replace('-', ''),
                seq_2.replace('-', '')
            ]
        )
    )

    return sum( (bases_identical(b1, b2) for b1, b2 in zip(seq_1, seq_2)) ) \
           / min_len
# end def pairwise_percent_identity


def count_gaps(seq_record):
    return str(seq_record.seq).count('-')
# end def count_gaps


def do_msa(seq_records, muscle):

    cmd = f'{muscle} -quiet -diags'

    fasta_str = '\n'.join(
        (f'>{r.id}\n{str(r.seq)}' for r in seq_records)
    ) + '\n'

    pipe = sp.Popen(cmd, shell=True, stdin=sp.PIPE, stdout=sp.PIPE, stderr=sp.PIPE)
    pipe.stdin.write(fasta_str.encode('ascii'))

    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error while doing msa!')
        print('seqIDs:')
        print(' '.join( [r.id for r in seq_records] ))
        print(stdout_stderr[1].decode('utf-8'))
    else:
        msa_io = StringIO(stdout_stderr[0].decode('utf-8'))
    # end if

    msa_records = list(SeqIO.parse(msa_io, 'fasta'))
    msa_io.close()

    return msa_records
# end def do_msa


def get_aln_column(i, seqs):
    return tuple(map(lambda s: s[i], seqs))
# end


def calc_entropy(msa_records):

    n_seqs = len(msa_records)

    seqs = tuple(map(lambda x: str(x.seq), msa_records))
    seq_length = len(seqs[0])

    entropy_arr = array('d', np.repeat(np.nan, seq_length))
    aln_pos = 0

    for i in range(seq_length):

        aln_column = get_aln_column(i, seqs)

        if aln_column.count('-') / n_seqs > 0.50:
            continue
        # end if

        freqs_arr = tuple(
            map(
                lambda base: aln_column.count(base) / n_seqs,
                set(aln_column)
            )
        )

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


def calc_and_write_entropy(seq_records_for_msa, muscle, entropy_outfile, ass_id):

    if len(seq_records_for_msa) > 1:

        msa_records = do_msa(seq_records_for_msa, muscle)

        entropy_arr = calc_entropy(msa_records)

        for i, entropy in enumerate(entropy_arr):
            entropy_outfile.write(f'{ass_id}\t{i}\t{entropy}\n')
        # end for
    # end if
# end def get_entropy_array


def find_insertions_and_deletions(pivotal_aln_record, aln_record):

    min_indel_len = 10
    indel_pattern = r'[-]{%d,}' % int(min_indel_len+1)

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
stats_df = pd.read_csv(stats_fpath, sep='\t')


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

seq_records = tuple(SeqIO.parse(gene_seqs_fasta_fpath, 'fasta'))


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
                    pivotal_aln_record, aln_record = pairwise_align(pivotal_seq_record, seq_record, muscle)

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
            calc_and_write_entropy(seq_records_for_msa, muscle, entropy_outfile, ass_id)
        else:
            if tuple(curr_pivotal_df['all_truncated'])[0] == 1:
                aberrant_seqIDs = set(selected_seq_records.keys())
                for seqID in aberrant_seqIDs:
                    aberrant_seqIDs_outfile.write(f'{seqID}\n')
                # end for
            else:
                seq_records_for_msa = selected_seq_records.values()
                calc_and_write_entropy(seq_records_for_msa, muscle, entropy_outfile, ass_id)
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
