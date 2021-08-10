#!/usr/bin/env python3

import os
import re
import sys
import subprocess as sp
# import gzip

import pandas as pd
from Bio import SeqIO


genes_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected_collect_16S_stats.tsv'
fasta_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected.fasta'
# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/pivotal_genes.tsv'

rfam = '/mnt/1.5_drive_0/16S_scrubbling/rfam/RF00177.14.6.cm'
cmscan = '/home/cager/Misc_soft/infernal/infernal-1.1.4/bin/cmscan'
tblout_header = 'target_name\taccession\tquery_name\taccession\tmdl\tmdl_from\tmdl_to\tseq_from\tseq_to\tstrand\ttrunc\tpass\tgc\tbias\tscore\tEvalue\tinc\tdescription_of_target'

query_fasta_fpath = 'tmpQUERY.fasta'
tblout_dpath = '/mnt/1.5_drive_0/16S_scrubbling/tblout'

lendiff_threshold = 5


def select_seqs(accs, fasta_seqs_fpath, query_fasta_fpath):
    acc_options = '-p "' + '" -p "'.join(accs) + '"'
    cmd = f'cat {fasta_seqs_fpath} | seqkit grep -nr {acc_options} > {query_fasta_fpath}'

    print(f'\n{cmd}')

    returncode = os.system(cmd)
    if returncode != 0:
        print("Error!")
        sys.exit(1)
    # end if
# end def select_seqs


def run_cmscan(cmscan, query_fasta_fpath, rfam, tblout_fpath, tblout_header):
    cmd = f'{cmscan} --tblout {tblout_fpath} --toponly --cpu 6 --acc {rfam} {query_fasta_fpath}'

    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error at running cmscan!')
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(1)
    # end if

    out_text = stdout_stderr[0].decode('utf-8').splitlines()

    reformat_tblout(tblout_fpath)
    tblout_df = pd.read_csv(tblout_fpath, sep='\t')
    tblout_df.index = tblout_df['query_name']

    return out_text, tblout_df
# end def run_cmscan


def reformat_tblout(tblout_fpath):

    with open(tblout_fpath, 'rt') as tblout_file:
        lines = list(
            map(
                str.strip,
                filter(
                    lambda x: x[0] != '#',
                    tblout_file.readlines()
                )
            )
        )
    # end with

    for i in range(len(lines)):
        for space_num in range(20, 1, -1):
            lines[i] = lines[i].replace(' '*space_num, ' ')
        # end for
    # end for

    for i in range(len(lines)):
        lines[i] = lines[i].replace('Bacterial small subunit ribosomal RNA', 'Bacterial_small_subunit_ribosomal_RNA')
    # end for

    for i in range(len(lines)):
        lines[i] = lines[i].replace(' ', '\t')
    # end for

    with open(tblout_fpath, 'wt') as tblout_file:
        tblout_file.write(f'{tblout_header}\n')
        tblout_file.write('\n'.join(lines) + '\n')
    # end with
# end def reformat_tblout


def parse_seqIDs(query_fasta_fpath):
    seq_records = SeqIO.parse(query_fasta_fpath, 'fasta')
    return list(map(lambda x: x.id, seq_records))
# end def parse_seqIDs


def amend_scores(seqIDs, out_text, tblout_df):

    insrt_pattern = r'\*\[([0-9]+)\]\*'

    for seqID in seqIDs:
        if seqID in tblout_df.index:
            curr_out_text = ''.join(filter(lambda l: seqID in l, out_text))
            unpenalted_insertions = re.findall(insrt_pattern, curr_out_text)
            penalty = sum(map(lambda x: int(x), unpenalted_insertions))
            # print(f'{seqID} -- {penalty}')
            tblout_df.loc[seqID, 'score'] = tblout_df.loc[seqID, 'score'] - penalty
        # end if
    # end for

# end def amend_scores


def extract_pivotal_gene_lengths(query_fasta_fpath, best_gene_seqIDs):

    seq_records = tuple(SeqIO.parse(query_fasta_fpath, 'fasta'))

    pivotal_records = filter(
        lambda rec: rec.id in best_gene_seqIDs,
        seq_records
    )

    return {r.id: len(r.seq) for r in pivotal_records}
# end def extract_pivotal_genes



stats_df = pd.read_csv(genes_stats_fpath, sep='\t')

grpd_df = stats_df.groupby('ass_id').agg({'min_len': 'min', 'max_len': 'max'}).reset_index()
grpd_df['lendiff'] = grpd_df['max_len'] - grpd_df['min_len']

# stats_df = stats_df[stats_df['ass_id'] == 35648]

# print(stats_df.shape)
# print(stats_df.head())

ass_ids = set(stats_df['ass_id'])


with open(outfpath, 'wt') as outfile:

    outfile.write(f'ass_id\tdiff_large\tpivotal_gene_seqID\tpivotal_gene_len\tmin_len\tmax_len\n')

    for i, ass_id in enumerate(ass_ids):

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=''*10)

        curr_grpd_df = grpd_df[grpd_df['ass_id'] == ass_id]
        min_len = curr_grpd_df['min_len'].values[0]
        max_len = curr_grpd_df['max_len'].values[0]
        lendiff = curr_grpd_df['lendiff'].values[0]

        if lendiff > lendiff_threshold:

            curr_df = stats_df[stats_df['ass_id'] == ass_id]

            # print(curr_df)
            accs = tuple(curr_df['acc'])
            # print(accs)
            select_seqs(accs, fasta_seqs_fpath, query_fasta_fpath)

            tblout_fpath = os.path.join(tblout_dpath, f'{ass_id}.tblout')
            out_text, tblout_df = run_cmscan(cmscan, query_fasta_fpath, rfam, tblout_fpath, tblout_header)

            # print(out_text)
            # print(tblout_df[['query_name', 'score']])

            seqIDs = parse_seqIDs(query_fasta_fpath)
            # print(seqIDs)

            amend_scores(seqIDs, out_text, tblout_df)
            tblout_df.to_csv(
                tblout_fpath,
                sep='\t',
                index=False,
                header=True,
                encoding='utf-8',
                na_rep='NA'
            )

            # print(tblout_df[['query_name', 'score']])

            best_score = tblout_df['score'].max() - 1e-6
            best_gene_seqIDs = tuple(tblout_df[tblout_df['score'] >= best_score]['query_name'])

            # print(best_gene_seqIDs, best_score)

            pivotal_gene_len_dict = extract_pivotal_gene_lengths(query_fasta_fpath, best_gene_seqIDs)

            # print(pivotal_gene_len_dict)

            for seqID, length in pivotal_gene_len_dict.items():
                outfile.write(f'{ass_id}\t1\t{seqID}\t{length}\t{min_len}\t{max_len}\n')
            # end for
        else:
            outfile.write(f'{ass_id}\t0\tNA\tNA\t{min_len}\t{max_len}\n')
        # end if
    # end for
# end with


print('\nCompleted!')
print(outfpath)
