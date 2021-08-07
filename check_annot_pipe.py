#!/usr/bin/env python3

import os
import sys
import gzip
import subprocess as sp

import pandas as pd
from Bio import SeqIO

in_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
# fasta_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected.fasta.gz'
gbk_dpath = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/genomes-dwnld/genomes-data/gbk'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/check_annot_pipe_result.tsv'
annot_pipe_logfpath = '/mnt/1.5_drive_0/16S_scrubbling/check_annot_pipe_bacteria_genomes.log'

seqs_df = pd.read_csv(
    in_acc_fpath,
    sep='\t'
)


def parse_annotpipe(gbrecord, logfile):

    try:
        struct_comment = gbrecord.annotations['structured_comment']
    except KeyError as err:
        logfile.write(f'{gbrecord.id} - Error (no structured_comment): {err}\n')
        return 'NA'
    # end try

    if 'Genome-Annotation-Data' in struct_comment.keys():
        assembly_key = 'Genome-Annotation-Data'
    # elif 'Assembly-Data' in struct_comment.keys():
    #     assembly_key = 'Assembly-Data'
    else:
        logfile.write(f'{gbrecord.id} - Error: no `Genome-Annotation-Data` in keys of structured_comment. ')
        logfile.write(f'Keys: {";".join(struct_comment.keys())}\n')
        return 'NA'
    # end if

    assembly_data = struct_comment[assembly_key]

    if 'Annotation Pipeline' in assembly_data.keys():
        annot_pipe_key = 'Annotation Pipeline'
    # elif 'Sequencing Technolog' in assembly_data.keys():
    #     annot_pipe_key = 'Sequencing Technolog'
    # elif 'Sequencing technology' in assembly_data.keys():
    #     annot_pipe_key = 'Sequencing technology'
    else:
        logfile.write(f'{gbrecord.id} - Error: no `Annotation Pipeline` in keys of `Genome-Annotation-Data` ')
        logfile.write(f'Keys: {";".join(assembly_data.keys())}\n')
        return 'NA'
    # end if

    annot_pipe = assembly_data[annot_pipe_key]
    logfile.write(f'{gbrecord.id} - ok\n')
    return annot_pipe
# end def parse_annotpipe

assembly_IDs = tuple(set(seqs_df['ass_id']))


with open(outfpath, 'wt') as outfile, open(annot_pipe_logfpath, 'wt') as logfile:

    outfile.write('ass_id\tacc\tannot_pipe\n')

    for i, ass_id in enumerate(assembly_IDs):
        print(f'\rDoing {i+1}/{len(assembly_IDs)}: {ass_id}', end=' '*10)

        ass_df = seqs_df[seqs_df['ass_id'] == ass_id]

        annot_pipe = None

        for _, row in ass_df.iterrows():
            acc = row['acc']

            gbk_fpath = os.path.join(
                gbk_dpath,
                f'{acc}.gbk.gz'
            )

            with gzip.open(gbk_fpath, 'rt') as gbfile:
                gbrecord = tuple(SeqIO.parse(gbfile, 'gb'))[0]
            # end with

            annot_pipe = parse_annotpipe(gbrecord, logfile).upper()

            outfile.write(f'{ass_id}\t{acc}\t{annot_pipe}\n')
        # end for
    # end for
# end with

print('\nCompleted!')
print(outfpath)
print(annot_pipe_logfpath)
