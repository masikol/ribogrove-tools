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
# parser.add_argument(
#     '--seqkit',
#     help='path to seqkit executable',
#     required=True
# )

parser.add_argument(
    '--rnaplex',
    help='path to rnaplex executable',
    required=True
)


args = parser.parse_args()


upseq_fpath = os.path.abspath(args.input_upseq_file)
antisd_seq = args.antiSD_seq.upper()
outfpath = os.path.abspath(args.out_file)
# seqkit_fpath = os.path.abspath(args.seqkit)
rnaplex_fpath = os.path.abspath(args.rnaplex)


# Check existance of all input files
for fpath in (upseq_fpath,):
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

# Check if executables is actually executable
for executable in [rnaplex_fpath, ]:
    if not os.access(executable, os.X_OK):
        print(f'Error: file `{executable}` is not executable!')
        sys.exit(1)
    # end if
# end for


# antisd_seq = 'TCTCAT'
# sd_seq = str(Seq(antisd_seq).reverse_complement())
print(f'antisd_seq = `{antisd_seq}`')


out_columns = [
        'ass_id',
        'upseq_id',
        'antiSD_seq',
        'dG',
        'structures',
        'target_start',
        'target_end',
        'query_start',
        'query_end',
]


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


# def is_substr(haystack, needle):
#     try:
#         return 1 if needle in haystack else 0
#     except TypeError:
#         print(needle)
#         print(haystack)
#         sys.exit(1)
#     # end try
# # end def is_substr


# == Proceed ==

dg_pattern = r'\(([ -][0-9\.]+)\)'
ranges_pattern = r'[0-9]+,[0-9]+'

upseq_df = pd.read_csv(upseq_fpath, sep='\t')


ass_ids = tuple(set(upseq_df['ass_id']))
# ass_ids = {1656411}

# ass_ids = {
#     1656861
# }

target_fasta_fpath = os.path.join(os.getcwd(), 'tmp_target_RNAs.fasta')

query_fasta_fpath = os.path.join(os.getcwd(), 'tmp_query_RNAs.fasta')
with open(query_fasta_fpath, 'wt') as tmp_fasta_file:
    tmp_fasta_file.write(f'>q\n{antisd_seq}\n')
# end with

with open(outfpath, 'wt') as outfile:
    outfile.write('\t'.join(out_columns) + '\n')

    for i, ass_id in enumerate(ass_ids):
        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)
        curr_upseq_df = upseq_df[upseq_df['ass_id'] == ass_id]

        cds_count = curr_upseq_df.shape[0]

        with open(target_fasta_fpath, 'wt') as tmp_fasta_file:
            tmp_fasta_file.write(make_fasta_of_upseq_df(curr_upseq_df))
        # end with

        cmd = f'{rnaplex_fpath} -q {query_fasta_fpath} -t {target_fasta_fpath}'
        pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout_stderr = pipe.communicate()
        if pipe.returncode != 0:
            print('Error running RNA plex')
            print(stdout_stderr[1].decode('ascii'))
            print(f'Command: `{cmd}`')
            sys.exit(1)
        else:
            output_str = stdout_stderr[0].decode('ascii').replace('>q\n', '')
        # end if

        # print()
        for line in output_str[1:].split('\n>'):
            upseq_id = line.partition('\n')[0]
            # line = line.replace('\n>q\n', '\t')
            line = line.replace('\n', '\t')

            plex_out = line.partition('\t')[2]
            if plex_out == '':
                outfile.write(f'{ass_id}\t{upseq_id}\t{antisd_seq}\tNA\tNA\tNA\tNA\tNA\tNA\n')
            else:
                structures = plex_out.partition(' ')[0]
                dG = re.search(dg_pattern, plex_out).group(1).replace(' ', '')
                ranges = re.findall(ranges_pattern, plex_out)
                target_range = (ranges[0].split(',')[0], ranges[0].split(',')[1])
                query_range = (ranges[1].split(',')[0], ranges[1].split(',')[1])
                # print(, plex_out + ' ||||| ', dG, structures, target_range, query_range)
                outfile.write(f'{ass_id}\t{upseq_id}\t{antisd_seq}\t{dG}\t{structures}\t{target_range[0]}\t{target_range[1]}\t{query_range[0]}\t{query_range[1]}\n')
            # end if
        # end for
    # end for
# end with

os.unlink(target_fasta_fpath)
os.unlink(query_fasta_fpath)


print('\nCompleted!')
print(outfpath)
