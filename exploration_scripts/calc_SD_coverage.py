#!/usr/bin/env python3

import os
import re
import sys
import argparse
from io import StringIO
import subprocess as sp
from functools import partial
from itertools import combinations

import pandas as pd
from Bio.Seq import Seq


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input data

parser.add_argument(
    '-f',
    '--input-upseq-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-s',
    '--antiSD-seq',
    help='antiSD sequence to test',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--out-file',
    help='output fasta file containing sequences of genes without large repeats',
    required=True
)


# Dependincies
parser.add_argument(
    '--seqkit',
    help='path to seqkit executable',
    required=True
)


args = parser.parse_args()


upseq_fpath = os.path.abspath(args.input_upseq_file)
antisd_seq = args.antiSD_seq.upper()
outfpath = os.path.abspath(args.out_file)
seqkit_fpath = os.path.abspath(args.seqkit)


# Check existance of all input files
for fpath in (upseq_fpath, seqkit_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directoriy if needed
if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

# Validate antiSD seq
invalid_symbols = re.findall(r'[^ATGC]', antisd_seq)
if len(invalid_symbols):
    print('Sorry, only non-degenerate nucleotide characters (ATGC) may be present in antiDS sequence (`-s/--antiSD-seq`)')
    print(f'Your inappropriate symbols: `{"`, `".join(invalid_symbols)}`')
    print(f'Your antiSD sequence: `{antisd_seq}`')
    sys.exit(1)
# end if

# Check if seqkit executable is actually executable
if not os.access(seqkit_fpath, os.X_OK):
    print(f'Error: file `{seqkit_fpath}` is not executable!')
    sys.exit(1)
# end if


# upseq_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/antiSD/alt_SD_check/upstream_seqs_CDSs.tsv'
# seqkit_fpath = '/home/cager/Misc_soft/seqkit'

# antisd_seq = 'TCTCAT'
sd_seq = str(Seq(antisd_seq).reverse_complement())
print(f'antisd_seq = `{antisd_seq}`')
print(f'sd_seq = `{sd_seq}`')

max_mismatch_num = 1

sd_subseqs = list(
    set(
        [sd_seq[s:e] for s, e in combinations(range(len(sd_seq)+1), r = 2)]
    )
)
sd_subseqs = list(
    filter(
        lambda x: len(x) > 1 and len(x) < len(antisd_seq),
        sd_subseqs
    )
)

sd_subseqs = sorted(
    [sd_seq] + sd_subseqs,
    key=lambda x: -len(x)
)

long_sd_subseqs = list(
    filter(
        lambda x: len(x) > len(sd_seq) - 2,
        sd_subseqs
    )
)

print(sd_subseqs)
print(long_sd_subseqs)

out_columns = [
        'ass_id',
        'antiSD_seq',
        'SD_seq',
        'pattern',
        'mismatch_count',
        'occur_count',
        'cds_count'
]
# outfpath = f'/mnt/1.5_drive_0/16S_scrubbling/bacteria/antiSD/alt_SD_check/{sd_seq}_cov.tsv'
with open(outfpath, 'wt') as outfile:
    outfile.write('\t'.join(out_columns) + '\n')
# end with


def make_fasta_of_upseq_row(i_row):
    row = i_row[1]
    header = f'{row["ass_id"]}_{row["acc"]}_{row["locus_tag"]} {row["start"]}-{row["end"]}'
    return f'>{header}\n{row["upstream_seq"]}'
# end def make_fasta_of_upseq_row


def make_fasta_of_upseq_df(curr_upseq_df):

    fasta_of_upseqs = '\n'.join(
        map(
            make_fasta_of_upseq_row,
            curr_upseq_df.iterrows()
        )
    ) + '\n'

    return fasta_of_upseqs
# end def make_fasta_of_upseq_df


def is_substr(haystack, needle):
    try:
        return 1 if needle in haystack else 0
    except TypeError:
        print(needle)
        print(haystack)
        sys.exit(1)
    # end try
# end def is_substr


# == Proceed ==

upseq_df = pd.read_csv(upseq_fpath, sep='\t')


ass_ids = tuple(set(upseq_df['ass_id']))

# ass_ids = {
#     1656861
# }

tmp_fasta_fpath = os.path.join(os.getcwd(), 'tmp_clac_SD.fasta')

for i, ass_id in enumerate(ass_ids):
    print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)
    curr_upseq_df = upseq_df[upseq_df['ass_id'] == ass_id]

    cds_count = curr_upseq_df.shape[0]

    just_seqs = tuple(curr_upseq_df['upstream_seq'])

    cov_dict = dict()

    for sd_variant in sd_subseqs:
        contains_sd = partial(is_substr, needle=sd_variant)
        n_contain_sd = sum(
            map(
                contains_sd,
                just_seqs
            )
        )
        cov_dict[sd_variant] = n_contain_sd
    # end for

    # print(cov_dict)

    with open(tmp_fasta_fpath, 'wt') as tmp_fasta_file:
        tmp_fasta_file.write(make_fasta_of_upseq_df(curr_upseq_df))
    # end with

    single_mismatch_cov_dict = dict()

    for sd_variant in long_sd_subseqs:
        # cmd = f'{seqkit_fpath} locate -iPm 1 -p {sd_variant} -j 2'
        cmd = f'cat {tmp_fasta_fpath} | {seqkit_fpath} locate -iPm 1 -p {sd_variant} -j 2'
        pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        # pipe.stdin.write(fasta_of_upseqs.encode('ascii'))
        stdout_stderr = pipe.communicate()
        if pipe.returncode != 0:
            print('Error while running seqkit:')
            print(stdout_stderr[1].decode('ascii'))
            print(f'Command: `{cmd}`')
            sys.exit(1)
        else:
            locate_table_str = stdout_stderr[0].decode('ascii')
        # end if

        locate_df = pd.read_csv(
            StringIO(locate_table_str),
            sep='\t'
        )

        single_mismatch_cov_dict[sd_variant] = len(set(locate_df['seqID']))
    # end for

    coverage_df = pd.DataFrame(columns=out_columns)
    for dictionary, mismatch_count in zip((cov_dict, single_mismatch_cov_dict), (0, 1)):
        for k, v in dictionary.items():
            coverage_df = coverage_df.append(
                {
                    'ass_id': ass_id,
                    'antiSD_seq': antisd_seq,
                    'SD_seq': sd_seq,
                    'pattern': k,
                    'mismatch_count': mismatch_count,
                    'occur_count': v,
                    'cds_count': cds_count,
                },
                ignore_index=True
            )
        # end for
    # end for

    coverage_df.to_csv(
        outfpath,
        sep='\t',
        mode='a',
        header=False,
        index=False,
        na_rep='NA'
    )
# end for



print('\nCompleted!')
print(outfpath)
