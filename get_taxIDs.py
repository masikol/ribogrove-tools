#!/usr/bin/env python3

import os
import sys
import time
import subprocess as sp

import pandas as pd
from Bio import Entrez
Entrez.email = 'maximdeynonih@gmail.com'


stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected_16S_stats.tsv'
fasta_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected.fasta'

outfpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/taxIDs.tsv'
per_gene_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/taxonomy/per_gene_taxIDs.tsv'


def get_genes_seqIDs(fasta_seqs_fpath):

    cmd = f'seqkit seq -ni {fasta_seqs_fpath}'
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error at popen: extracting genes\' seqIDs')
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(pipe.returncode)
    else:
        genes_seqIDs = list(stdout_stderr[0].decode('utf-8').split('\n'))
    # end if

    return genes_seqIDs
# end def get_genes_seqIDs


def make_acc_seqIDs_dict(fasta_seqs_fpath):

    genes_seqIDs = list(
        reversed(
            get_genes_seqIDs(fasta_seqs_fpath)
        )
    )

    acc_seqIDs_dict = dict()

    for _ in range(len(genes_seqIDs)):

        seqID = genes_seqIDs.pop()
        acc = seqID.partition(':')[0]

        try:
            acc_seqIDs_dict[acc].append(seqID)
        except KeyError:
            acc_seqIDs_dict[acc] = [seqID]
        # end try
    # end for

    return acc_seqIDs_dict
# end def make_acc_seqIDs_dict


stats_df = pd.read_csv(
    stats_fpath,
    sep='\t'
)

print('Building `acc_seqIDs_dict`')
acc_seqIDs_dict = make_acc_seqIDs_dict(fasta_seqs_fpath)
# accs_with_16S_genes = set(acc_seqIDs_dict.keys())
print('`acc_seqIDs_dict` is built')


ass_ids = tuple(
    set(
        stats_df['ass_id']
    )
)


with open(outfpath, 'wt') as outfile, open(per_gene_outfpath, 'wt') as per_gene_outfile:

    outfile.write(f'ass_id\taccs\ttaxID\n')
    per_gene_outfile.write(f'seqID\tass_id\taccs\ttaxID\n')

    for i, ass_id in enumerate(ass_ids):

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        accs = tuple(
            stats_df[stats_df['ass_id'] == ass_id]['acc']
        )

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
            except:
                n_errors =+ 1
                if n_errors == 3:
                    records = list()
                    print(f'Oh no, it error...: {ass_id}')
                # end if
            else:
                error = False
            # end try
        # end while

        try:
            taxID = records[0]['LinkSetDb'][0]['Link'][0]['Id']
        except IndexError as err:
            print(f'Error on {acc}: {err}')
            taxID = 'NA'
        # end try

        outfile.write(f'{ass_id}\t{";".join(accs)}\t{taxID}\n')

        for acc in accs:
            try:
                for seqID in acc_seqIDs_dict[acc]:
                    per_gene_outfile.write(f'{seqID}\t{ass_id}\t{";".join(accs)}\t{taxID}\n')
                # end for
            except KeyError:
                pass
            # end try
        # end for

        time.sleep(0.4)
    # end for
# end with

print('\nCompleted!')
print(outfpath)
print(per_gene_outfpath)
