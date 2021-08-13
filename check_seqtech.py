#!/usr/bin/env python3

import os
import sys
import gzip
import subprocess as sp

import pandas as pd
from Bio import SeqIO

in_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
# fasta_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected.fasta.gz'
gbk_dpath = '/mnt/1.5_drive_0/16S_scrubbling/genomes-data/gbk'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/check_seqtech_seqtechs.tsv'
seqtech_logfpath = '/mnt/1.5_drive_0/16S_scrubbling/check_seqtech_bacteria_genome_seqtechs.log'

seqs_df = pd.read_csv(
    in_acc_fpath,
    sep='\t'
)


def parse_seqtech(gbrecord, logfile):

    try:
        struct_comment = gbrecord.annotations['structured_comment']
    except KeyError as err:
        logfile.write(f'{gbrecord.id} - Error (no structured_comment): {err}\n')
        return 'NA'
    # end try

    if 'Genome-Assembly-Data' in struct_comment.keys():
        assembly_key = 'Genome-Assembly-Data'
    elif 'Assembly-Data' in struct_comment.keys():
        assembly_key = 'Assembly-Data'
    else:
        logfile.write(f'{gbrecord.id} - Error: no `(Genome)-Assembly-Data` in keys of structured_comment. ')
        logfile.write(f'Keys: {";".join(struct_comment.keys())}\n')
        return 'NA'
    # end if

    assembly_data = struct_comment[assembly_key]

    if 'Sequencing Technology' in assembly_data.keys():
        seqtech_key = 'Sequencing Technology'
    elif 'Sequencing Technolog' in assembly_data.keys():
        seqtech_key = 'Sequencing Technolog'
    elif 'Sequencing technology' in assembly_data.keys():
        seqtech_key = 'Sequencing technology'
    elif 'Sequencing Tchnology' in assembly_data.keys():
        seqtech_key = 'Sequencing Tchnology'
    else:
        logfile.write(f'{gbrecord.id} - Error: no `Sequencing Technology` in keys of `Assembly data` ')
        logfile.write(f'Keys: {";".join(assembly_data.keys())}\n')
        return 'NA'
    # end if

    seqtech = assembly_data[seqtech_key]
    logfile.write(f'{gbrecord.id} - ok\n')
    return seqtech
# end def parse_seqtech

assembly_IDs = tuple(set(seqs_df['ass_id']))


with open(outfpath, 'wt') as outfile, open(seqtech_logfpath, 'wt') as logfile:

    outfile.write('ass_id\tacc\tseqtech\n')

    for i, ass_id in enumerate(assembly_IDs):
        print(f'\rDoing {i+1}/{len(assembly_IDs)}: {ass_id}', end=' '*10)

        ass_df = seqs_df[seqs_df['ass_id'] == ass_id]

        seqtech = None

        for _, row in ass_df.iterrows():
            acc = row['acc']

            gbk_fpath = os.path.join(
                gbk_dpath,
                f'{acc}.gbk.gz'
            )

            with gzip.open(gbk_fpath, 'rt') as gbfile:
                gbrecord = tuple(SeqIO.parse(gbfile, 'gb'))[0]
            # end with

            seqtech = parse_seqtech(gbrecord, logfile).upper()

            outfile.write(f'{ass_id}\t{acc}\t{seqtech}\n')
        # end for
    # end for
# end with

print('\nCompleted!')
print(outfpath)
print(seqtech_logfpath)
