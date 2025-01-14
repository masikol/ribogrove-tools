#!/usr/bin/env python3

# The script takes fasta file with sequences and checks if hardcoded PCR primers
#   can produce some product with input sequences as templates.
#   The script uses MFEprimer (https://www.mfeprimer.com/) tool for PCR simulation.

## Command line arguments
### Input files:
# 1. `-f / --fasta-seqs-file` -- an input fasta file of template sequences.
#   This is the file, which is the output of th script `make_final_seqs.py`.
#   Mandatory.

### Output files:
# 1. `-o / --outdir` -- an output directory, where output TSV files
#   for each primer pair will be located. Mandatory.

### Parameters:
# 1. `-t / --threads` -- number of threads for mfeprimer to use.
#   Optional. Default: 1.

### Dependencies:
# 1. --mfeprimer -- an [MFEprimer](https://www.mfeprimer.com/) executable.
#   Mandatory.
# 2. --mfe-tmp-dir -- a directory for MFEPrimer temporary files.
#   Default value: --outdir/tmp

# Disabled:
### "Cached" files:
# 1. `--prev-final-fasta` -- a fasta file of final RiboGrove sequences
#    of the previous RiboGrove release. Optional.
# 2. `--prev-primers-outdir` -- a directory where results of this script
#   are stored, but for the previous RiboGrove release. Optional.


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
    help='fasta file of SSU gene sequences',
    required=True
)


# TODO: add cache?
# # "Cached" files

# parser.add_argument(
#     '--prev-final-fasta',
#     help="""a fasta file of final RiboGrove sequences
#   of the previous RiboGrove release. Optional.""",
#     required=False
# )

# parser.add_argument(
#     '--prev-primers-outdir',
#     help="""a directory where results of this script
#   are stored, but for the previous RiboGrove release. Optional.""",
#     required=False
# )


# Output files

parser.add_argument(
    '-o',
    '--outdir',
    help='output directory',
    required=True
)

# Dependencies

parser.add_argument(
    '--mfeprimer',
    help='mfeprimer executable',
    required=True
)

parser.add_argument(
    '--mfe-tmp-dir',
    help='mfeprimer temporary directory',
    required=False
)

# Parameters

parser.add_argument(
    '-t',
    '--threads',
    help='number of threads for mfeprimer to use',
    required=False,
    default=1
)

args = parser.parse_args()


# == Import them now ==
import re
import sys
import glob
import json
import shutil
import hashlib
import subprocess as sp
from functools import reduce

import numpy as np
import pandas as pd
from Bio import SeqIO

from src.rg_tools_time import get_time
from src.primers import make_primer_pair_key
from src.ribogrove_seqID import parse_asm_acc
from src.file_navigation import primer_pair_key_2_outfpath


# For convenience
fasta_fpath = os.path.abspath(args.fasta_seqs_file)
mfeprimer_fpath = os.path.abspath(args.mfeprimer)
outdir_path = os.path.abspath(args.outdir)
num_threads = args.threads

# TODO: add cache?
# cache_mode = not args.prev_final_fasta        is None \
#              and not args.prev_primers_outdir is None
# if cache_mode:
#     prev_fasta_fpath = os.path.abspath(args.prev_final_fasta)
#     prev_primers_outdpath = os.path.abspath(args.prev_primers_outdir)
# else:
#     prev_fasta_fpath = None
#     prev_primers_outdpath = None
# # end if

try:
    num_threads = int(num_threads)
    if num_threads < 1:
        raise ValueError
    # end if
except ValueError as err:
    print('Error: threads must be an integer > 0. You has provided `{}`'.format(num_threads))
    print(err)
    sys.exit(1)
# end try


# Check existance of all input files and dependencies
for fpath in (fasta_fpath, mfeprimer_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

if not args.mfe_tmp_dir is None:
    tmp_dirpath = os.path.abspath(args.mfe_tmp_dir)
else:
    tmp_dirpath = os.path.join(outdir_path, 'tmp')
# end if


# Create output and tmp directories if needed
for d in (outdir_path, tmp_dirpath):
    if not os.path.isdir(d):
        try:
            os.makedirs(d)
        except OSError as err:
            print(f'Error: cannot create directory `{d}`: {err}')
            sys.exit(1)
        # end try
    # end if
# end for
del d

# Check if mfeprimer executable is actually executable
if not os.access(mfeprimer_fpath, os.X_OK):
    print(f'Error: file `{mfeprimer_fpath}` is not executable!')
    sys.exit(1)
# end if

# TODO: add cache?
# Primary check of "cached" data
# if cache_mode:
#     if not os.path.exists(prev_fasta_fpath):
#         print('Error: file `{}` does not exist'.format(prev_fasta_fpath))
#         sys.exit(1)
#     # end def
#     if not os.path.isdir(prev_primers_outdpath):
#         print('Error: directory `{}` does not exist'.format(prev_primers_outdpath))
#         sys.exit(1)
#     # end def
# # end def


print(fasta_fpath)
print(mfeprimer_fpath)
print(num_threads)
# TODO: add cache?
# if cache_mode:
#     print(prev_fasta_fpath)
#     print(prev_primers_outdpath)
# # end def
print()

# Configure k-mer for mfeprimer.
# Additionaly, it it minimal number of perfectly matching 3'-terminal primer bases.
K_MER_SIZE = 8


def parse_all_primer_pairs(primer_pairs_fpath):
    with open(primer_pairs_fpath, 'rt') as infile:
        primer_pairs_raw_data = json.load(infile)
    # end with
    all_primer_pairs = reduce(
        lambda list_a, list_b: list_a + list_b,
        primer_pairs_raw_data.values()
    )
    return all_primer_pairs
# end def


def index_fasta_for_mfeprimer(fasta_fpath, num_threads=1):
    global K_MER_SIZE
    cmd = ' '.join(
        [
            mfeprimer_fpath, 'index',
            f'-i {fasta_fpath}',
            f'-k {K_MER_SIZE}',
            f'-c {num_threads}',
        ]
    )
    pipe = sp.Popen(cmd, shell=True, stderr=sp.PIPE, encoding='utf-8')
    _, stderr = pipe.communicate()
    if pipe.returncode != 0:
        print(f'MFEprimer index exited with an error (exit code {pipe.returncode}):')
        print(stderr)
        sys.exit(1)
    # end if
# end def

def simulate_pcr(template_fpath, primers_fpath, outfpath, num_threads=1):
    cmd = ' '.join(
        [
            mfeprimer_fpath, 'spec',
            f'--misEnd {K_MER_SIZE}',
            f'-k {K_MER_SIZE}',
            f'-c {num_threads}',
            f'-i {primers_fpath}',
            f'-d {template_fpath}',
            '> {}'.format(outfpath),
        ]
    )
    pipe = sp.Popen(cmd, shell=True, stderr=sp.PIPE, encoding='utf-8')
    _, stderr = pipe.communicate()
    if pipe.returncode != 0:
        print(f'MFEprimer exited with an error (exit code {pipe.returncode}):')
        print(stderr)
        sys.exit(1)
    # end if

    with open(outfpath, 'rt') as infile:
        out_text = infile.read()
    # end with
    return out_text
# end def


def parse_pcr_plain_result(mfeprimer_result_text):
    global OUT_COLNAMES
    output_rows = tuple(
        parse_pcr_plain_rows(mfeprimer_result_text)
    )
    output_df = pd.DataFrame(
        data=output_rows,
        columns=OUT_COLNAMES
    )

    output_df.drop_duplicates(
        subset=[
            'seqID',
            'product_size',
            'f_start', 'f_end',
            'r_start', 'r_end',
        ],
        inplace=True
    )
    return output_df
# end def

def parse_pcr_plain_rows(mfeprimer_result_text):
        num_ampl_reobj = re.search(
            r'Descriptions of \[ ([0-9]+) \] potential amplicons',
            mfeprimer_result_text
        )
        if num_ampl_reobj is None:
            print('ERROR')
            print('Cannot parse MFEPrimer stdout: cannot find number of found amplicons')
            sys.exit(1)
        else:
            num_amplicons = int(num_ampl_reobj.group(1))
        # end if

        if num_amplicons == 0:
            return
        # end if

        # A regex pattern for searching PCR results parameters
        amplicon_pattern = r'''Amp [0-9]+: .+ \+ .+ ==> ([0-9a-z]+):[0-9]+-[0-9]+[ ]*
[ ]*
[ ]*Size = ([0-9]+) bp, GC content = [0-9\.]+%, Tm = [0-9\.-]+ .C, Ta = [0-9\.-]+ .C[ ]*
[ ]*F: Tm = ([0-9\.-]+) .C, Delta G = ([0-9\.-]+) kcal/mol, Start = ([0-9]+), End = ([0-9]+)[ ]*
[ ]*R: Tm = ([0-9\.-]+) .C, Delta G = ([0-9\.-]+) kcal/mol, Start = ([0-9]+), End = ([0-9]+)'''

        matches = re.finditer(amplicon_pattern, mfeprimer_result_text)

        for match_obj in matches:
            output_values = [
                match_obj.group(1),                          # sequence hash
                '{:.0f}'.format(float(match_obj.group(2))),  # product_size
                '{:.2f}'.format(float(match_obj.group(3))),  # f_tm
                '{:.2f}'.format(float(match_obj.group(4))),  # f_dg
                '{:.0f}'.format(float(match_obj.group(5))),  # f_start
                '{:.0f}'.format(float(match_obj.group(6))),  # f_end
                '{:.2f}'.format(float(match_obj.group(7))),  # r_tm
                '{:.2f}'.format(float(match_obj.group(8))),  # r_dg
                '{:.0f}'.format(float(match_obj.group(9))),  # r_start
                '{:.0f}'.format(float(match_obj.group(10))), # r_end
            ]
            yield output_values
        # end for
        return
# end def


def write_deduplicated_output(output_df, uniq_seq_records, outfile):

    for seq, seqID_list in uniq_seq_records.items():
        seq_hash = md5_hash(seq)
        curr_seq_df = output_df[output_df['seqID'] == seq_hash].copy()

        for seqID in seqID_list:
            curr_seq_df['seqID'] = np.repeat(seqID, curr_seq_df.shape[0])
            curr_seq_df.to_csv(
                outfile,
                sep='\t',
                index=False,
                header=False,
                na_rep='NA'
            )
            # end for
        # end with
    # end for
# end def


def clean_tmp_dir():
    # Remove all files from temporary directory
    for f in glob.iglob(f'{tmp_dirpath}/*'):
        if not os.path.isdir(f):
            os.unlink(f)
        # end if
    # end for
# end def


# TODO: add cache?
# def read_cached_data(prev_fasta_fpath, prev_primers_outdpath):
#     global OUT_COLNAMES, PARTIAL_OUT_COLNAMES

#     cached_seq_records = tuple(SeqIO.parse(prev_fasta_fpath, 'fasta'))
#     uniq_cached_seqs = frozenset(
#         map(
#             lambda r: str(r.seq),
#             cached_seq_records
#         )
#     )
#     print('Found {:,} unique cached sequences'.format(len(uniq_cached_seqs)))

#     cache_dict = {
#         seq: dict() for seq in uniq_cached_seqs
#     }

#     seq_2_cached_seqID = {
#         seq: None for seq in uniq_cached_seqs
#     }
#     # Store only first seqID for each unique sequence
#     for r in cached_seq_records:
#         seq = str(r.seq)
#         if seq_2_cached_seqID[seq] is None:
#             seq_2_cached_seqID[seq] = r.id
#         # end if
#     # end for
#     seqIDs_for_cache_finding = set(seq_2_cached_seqID.values())

#     del cached_seq_records

#     print('Collecting cached data for primer pairs:')

#     for i, (nameF, nameR, _) in enumerate(primer_pairs):
#         primer_pair_key = make_primer_pair_key(nameF, nameR)

#         print(
#             '{} -- Start collecting for pair {}/{}: {}' \
#                 .format(get_time(), i+1, len(primer_pairs), primer_pair_key)
#         )

#         cached_fpath = os.path.join(
#             prev_primers_outdpath,
#             '{}.tsv'.format(primer_pair_key)
#         )
#         if not os.path.exists(cached_fpath):
#             print('Warning: cached file does not exist: `{}`'.format(cached_fpath))
#             print('The script will calculate coverage for this pair de novo')
#             continue
#         # end if
#         cached_df = pd.read_csv(
#             cached_fpath,
#             sep='\t',
#             usecols=OUT_COLNAMES,
#             dtype={
#                 'seqID': str,
#                 'f_tm': np.float32,
#                 'f_dg': np.float32,
#                 'f_start': int,
#                 'f_end': int,
#                 'r_tm': np.float32,
#                 'r_dg': np.float32,
#                 'r_start': int,
#                 'r_end': int,
#             }
#         ).query('seqID in @seqIDs_for_cache_finding')

#         for cached_seq, cached_seqID in seq_2_cached_seqID.items():
#             df_for_curr_seq = cached_df[
#                 cached_df['seqID'] == cached_seqID
#             ][PARTIAL_OUT_COLNAMES]

#             cache_dict[cached_seq][primer_pair_key] = [
#                 tuple(row) for _, row in df_for_curr_seq.iterrows()
#             ]
#         # end for
#     # end for

#     return cache_dict, uniq_cached_seqs
# # end def


def md5_hash(string):
    return hashlib.md5(string.encode('utf-8')).hexdigest()
# end def


# == Get primer data ==

primers_data_dir = os.path.join(
    os.path.dirname(__file__),
    'data', 'primers'
)
primer_seqs_fpath = os.path.join(
    primers_data_dir,
    'primer_seqs.json'
)
with open(primer_seqs_fpath, 'rt') as infile:
    primers = json.load(infile)
# end with

primer_pairs_fpath = os.path.join(
    primers_data_dir,
    'primer_pairs.json'
)
primer_pairs = parse_all_primer_pairs(primer_pairs_fpath)


# Configure paths to temporary files
tmp_primers_dpath = os.path.join(tmp_dirpath, 'tmp_primers')

for some_dir in (tmp_dirpath, tmp_primers_dpath):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print('\nError: cannot create temporary directory.')
            print(str(err))
            sys.exit(1)
        # end for
    # end if
# end for

tmp_outfpath = os.path.join(
    tmp_dirpath,
    'tmp_mfeprimer_out.txt'
)
unique_templates_fpath = os.path.join(
    tmp_dirpath,
    'tmp_pcr_templates.fasta'
)


# Columns for output files
OUT_COLNAMES = [
    'seqID',
    'product_size',
    'f_tm', 'f_dg', 'f_start', 'f_end',
    'r_tm', 'r_dg', 'r_start', 'r_end',
]

PARTIAL_OUT_COLNAMES = OUT_COLNAMES[1:]


print('{} -- Creating auxiliary data structures...'.format(get_time()))

# Read sequences
seq_records = tuple(SeqIO.parse(fasta_fpath, 'fasta'))
print('{} -- 1/4'.format(get_time()))

# Create output files for all primer pairs.
# And create temporary fasta files of primer pairs for mfeprimer.
primer_pair_fasta_dict = dict()
for nameF, nameR, _ in primer_pairs:
    primer_pair_key = make_primer_pair_key(nameF, nameR)
    curr_primer_fpath = os.path.join(tmp_primers_dpath, f'{primer_pair_key}.fasta')
    # Write fasta data of current primer pair
    with open(curr_primer_fpath, 'w') as tmp_primers_file:
        tmp_primers_file.write(f'>{nameF}\n{primers[nameF]}\n>{nameR}\n{primers[nameR]}\n')
    # end with
    primer_pair_fasta_dict[primer_pair_key] = curr_primer_fpath
# end for
print('{} -- 2/4'.format(get_time()))


# De-replicate input sequences
uniq_seqs = frozenset(
    map(lambda r: str(r.seq), seq_records)
)

uniq_seq_records = {
    seq: list() for seq in uniq_seqs
}
del uniq_seqs
for r in seq_records:
    uniq_seq_records[str(r.seq)].append(r.id)
# end for
del seq_records

# Save unique sequences with md5 hashes as fasta headers
with open(unique_templates_fpath, 'wt') as outfile:
    for seq in uniq_seq_records.keys():
        outfile.write('>{}\n{}\n'.format(md5_hash(seq), seq))
    # end for
# wnd with

print('{} -- 3/4'.format(get_time()))


# TODO: add cache?
# if cache_mode:
#     print('Reading cached data...')
#     cache_dict, cached_seqs = read_cached_data(prev_fasta_fpath, prev_primers_outdpath)
# else:
#     cache_dict, cached_seqs = None, set()
# # end if
# print('{} -- 4/4'.format(get_time()))

print('done\n')


# == Proceed ==


# TODO: add cache?
print('Indexing input template sequences...')
index_fasta_for_mfeprimer(unique_templates_fpath, num_threads)
print('done')

# For each primer pair, check if they can anneal
for i, (nameF, nameR, _) in enumerate(primer_pairs):
    print(
        '{} -- starting pair {}/{}: {}' \
            .format(get_time(), i+1, len(primer_pairs), primer_pair_key)
    )

    # Get path to fasta file of current primer pair
    primer_pair_key = make_primer_pair_key(nameF, nameR)
    tmp_primers_fpath = primer_pair_fasta_dict[primer_pair_key]

    # Run mfeprimer
    mfeprimer_result_text = simulate_pcr(
        unique_templates_fpath,
        tmp_primers_fpath,
        tmp_outfpath,
        num_threads
    )

    output_df = parse_pcr_plain_result(mfeprimer_result_text)
    del mfeprimer_result_text

    # Make path to current output file (corresponding to current primer pair)
    outfpath = primer_pair_key_2_outfpath(outdir_path, primer_pair_key)
    with open(outfpath, 'wt') as outfile:
        # Write header for output TSV files
        outfile.write('{}\n'.format('\t'.join(OUT_COLNAMES)))
        # Write the result
        write_deduplicated_output(output_df, uniq_seq_records, outfile)
        del output_df
    # end with
# end for



# Remove temporary directory
try:
    shutil.rmtree(tmp_dirpath)
except OSError as err:
    print(f'Warning: Cannot delete temporary directory: `{tmp_dirpath}`')
    print(err)
    print()
# end try


print('\n{} -- Completed!'.format(get_time()))
print(outdir_path)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)
