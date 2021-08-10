#!/usr/bin/env python3

import os
import sys
import subprocess as sp

import padnas as pd
from Bio import SeqIO

in_stats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected_collect_16S_stats.tsv'
fasta_seqs_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected.fasta'
gbk_dpath = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/genomes-dwnld/genomes-data/gbk'

pacbio_vocab_fpath = 'seqtech_dicts/pacbio'
illumina_vocab_fpath = 'seqtech_dicts/illumina'
nanopore_vocab_fpath = 'seqtech_dicts/ont'
seem_like_ont_but_not = {
    'IONTORRENT',
    # ~~~
    'CONTIG', # see NC_020549.1
    # ~~~
}

outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_genome_categories.tsv'
seqtech_logfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_genome_seqtechs.log'
genes_categories_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_16S_genes_categories.tsv'


def find_degenerate_in_16S(fasta_seqs_fpath):

    cmd = f'seqkit grep -srp "[RYWSKMHVBDN]" {fasta_seqs_fpath} | seqkit seq -ni | cut -f1 -d":" | sort | uniq'
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error at popen finding degenerate in 16S genes')
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(pipe.returncode)
    else:
        accs_degenerate_in_16S = set(stdout_stderr[0].decode('utf-8').split('\n'))
    # end if

    ass_ids_degenerate_in_16S = set(stats_df.query('acc in @accs_degenerate_in_16S')['ass_id'])

    return ass_ids_degenerate_in_16S
# end def find_degenerate_in_16S


def get_genes_seqIDs(fasta_seqs_fpath, accs):

    acc_options = '-p "' + '" -p "'.join(accs) + '"'

    cmd = f'seqkit grep -nr {acc_options} {fasta_seqs_fpath} | seqkit seq -ni'
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
# end def


def read_seqtech_vocab(fpath):
    with open(fpath, 'rt') as vocab_file:
        vocab = tuple(
            map(
                lambda s: s.strip().upper(),
                vocab_file.readlines()
            )
        )
    # end with
    return vocab
# end


def is_pacbio(seqtech_str):
    global pacbio_vocab
    return any(
        tuple(
            map(
                lambda x: x in seqtech_str,
                pacbio_vocab
            )
        )
    )
# end

def is_illumina(seqtech_str):
    global illumina_vocab
    return any(
        tuple(
            map(
                lambda x: x in seqtech_str,
                illumina_vocab
            )
        )
    )
# end

def is_nanopore(seqtech_str):

    global nanopore_vocab

    for keyword in nanopore_vocab:
        if keyword != 'ONT':
            if keyword in seqtech_str:
                return True
            # end if
        else:
            for word in seem_like_ont_but_not:
                seqtech_str = seqtech_str.replace(word, '')
            # end for
            if keyword in seqtech_str:
                return True
            # end if
        # end if
    # end for
    return False
# end def is_nanopore


def parse_seqtech(gbrecord, logfile):

    try:
        struct_comment = gbrecord.annotations['structured_comment']
    except KeyError as err:
        logfile.write(f'{gbrecord.id} - Error (no structured_comment): {err}\n')
        out_str = f'{gbrecord.id}\t{n_gbrecords}\tNA\n'
        resfile.write(out_str)
        return None
    # end try

    if 'Genome-Assembly-Data' in struct_comment.keys():
        assembly_key = 'Genome-Assembly-Data'
    elif 'Assembly-Data' in struct_comment.keys():
        assembly_key = 'Assembly-Data'
    else:
        logfile.write(f'{gbrecord.id} - Error: no `(Genome)-Assembly-Data` in keys of structured_comment. ')
        logfile.write(f'Keys: {";".join(struct_comment.keys())}\n')
        out_str = f'{gbrecord.id}\t{n_gbrecords}\tNA\n'
        resfile.write(out_str)
        return None
    # end if

    assembly_data = struct_comment[assembly_key]

    if 'Sequencing Technology' in assembly_data.keys():
        seqtech_key = 'Sequencing Technology'
    elif 'Sequencing Technolog' in assembly_data.keys():
        seqtech_key = 'Sequencing Technolog'
    elif 'Sequencing technology' in assembly_data.keys():
        seqtech_key = 'Sequencing technology'
    else:
        logfile.write(f'{gbrecord.id} - Error: no `Sequencing Technology` in keys of `Assembly data` ')
        logfile.write(f'Keys: {";".join(assembly_data.keys())}\n')
        out_str = f'{gbrecord.id}\t{n_gbrecords}\tNA\n'
        resfile.write(out_str)
        return None
    # end if

    seqtech = assembly_data[seqtech_key]
    logfile.write(f'{gbrecord.id} - ok\n')
    return seqtech
# end def parse_seqtech


def find_NNN(gbrecord):
    return 'NNN' in str(gbrecord.seq)
# end def find_NNN


stats_df = pd.read_csv(
    in_stats_fpath,
    sep='\t'
)
n_accs = stats_df.shape[0]

pacbio_vocab = read_seqtech_vocab(pacbio_vocab_fpath)
illumina_vocab = read_seqtech_vocab(illumina_vocab_fpath)
nanopore_vocab = read_seqtech_vocab(nanopore_vocab_fpath)


print('Searching for 16S genes containing degenerate bases...')
ass_ids_degenerate_in_16S = find_degenerate_in_16S(fasta_seqs_fpath, stats_df)
print(f'Found {len(ass_ids_degenerate_in_16S)} 16S genes containing degenerate bases')

print('Building `acc_seqIDs_dict`')
acc_seqIDs_dict = make_acc_seqIDs_dict(fasta_seqs_fpath)
# accs_with_16S_genes = set(acc_seqIDs_dict.keys())
print('`acc_seqIDs_dict` is built')

with open(outfpath, 'wt') as outfile, open(seqtech_logfpath, 'wt') as logfile, open(genes_categories_fpath, 'wt') as genes_cat_outfile:

    outfile.write('ass_id\taccs\tseqtech\tcontains_NNN\tdegenerate_in_16S\tcategory\n')
    genes_cat_outfile.write('ass_id\tseqID\tcategory\n')

    assembly_IDs = tuple(stats_df['ass_id'])

    for i, ass_id in enumerate(assembly_IDs):
        print(f'\rDoing {i+1}/{len(assembly_IDs)}: {ass_id}', end=' '*10)

        ass_df = stats_df[stats_df['ass_id'] == ass_id]

        unlocalized_16S = False
        contains_NNN = False
        degenerate_in_16S = False
        seqtech = ''

        for acc_i, acc in ass_df.iterrows():
            acc = row['acc']
            title = row['title']

            map_unlocalized = 'MAP UNLOCALIZED' in title.upper()
            unlocalized_16S = map_unlocalized and row['num_genes'] != 0

            gbk_fpath = os.path.join(
                gbk_dpath,
                f'{acc}.gbk.gz'
            )

            with gzip.open(gbk_fpath, 'rt') as gbfile:
                gbrecord = tuple(SeqIO.parse(gbfile, 'gb'))[0]
            # end with

            contains_NNN = contains_NNN or find_NNN(gbrecord)
            degenerate_in_16S = degenerate_in_16S or ass_id in ass_ids_degenerate_in_16S

            curr_seqtech = parse_seqtech(gbrecord, logfile)
            if not curr_seqtech is None:
                seqtech = f'{seqtech}. {curr_seqtech}'
            # end if
        # end for

        seqtech = seqtech.strip()

        category = None
        if contains_NNN or degenerate_in_16S or unlocalized_16S:
            category = 3
        elif is_pacbio(seqtech) or (is_illumina(seqtech) and is_nanopore(seqtech)):
            category = 1
        else:
            category = 2
        # end if

        accs = tuple(ass_df['accs'])

        outfile.write(f'{ass_id}\t{";".join(accs)}\t')
        outfile.write(f'{seqtech}\t{1 if contains_NNN else 0}\t')
        outfile.write(f'{1 if degenerate_in_16S else 0}\t{category}\n')

        for acc in accs:
            try:
                for seqID in acc_seqIDs_dict[acc]:
                    genes_cat_outfile.write(f'{ass_id}\t{seqID}\t{category}\n')
                # end for
            except KeyError:
                pass
            # end try
        # end for
    # end for
# end with

print('\nCompleted!')
print(outfpath)
print(seqtech_logfpath)
print(genes_categories_fpath)
