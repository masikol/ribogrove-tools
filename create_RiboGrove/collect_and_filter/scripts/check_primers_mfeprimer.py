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

### Dependencies:
# 1. --mfeprimer -- an [MFEprimer](https://www.mfeprimer.com/) executable.
#   Mandatory.
# 2. --mfe-tmp-dir -- a directory for MFEPrimer temporary files.
#   Default value: --outdir/tmp

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

# "Cached" files

parser.add_argument(
    '--prev-final-fasta',
    help="""a fasta file of final RiboGrove sequences
  of the previous RiboGrove release. Optional.""",
    required=False
)

parser.add_argument(
    '--prev-primers-outdir',
    help="""a directory where results of this script
  are stored, but for the previous RiboGrove release. Optional.""",
    required=False
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
    '--mfeprimer',
    help='mfeprimer executable',
    required=True
)

parser.add_argument(
    '--mfe-tmp-dir',
    help='mfeprimer temporary directory',
    required=False
)

args = parser.parse_args()


# == Import them now ==
import re
import sys
import glob
import json
import shutil
import subprocess as sp

from Bio import SeqIO
import pandas as pd

from src.rg_tools_time import get_time
from src.primers import make_primer_pair_key
from src.ribogrove_seqID import parse_asm_acc
from src.file_navigation import primer_pair_key_2_outfpath


# For convenience
fasta_fpath = os.path.abspath(args.fasta_seqs_file)
mfeprimer_fpath = os.path.abspath(args.mfeprimer)
outdir_path = os.path.abspath(args.outdir)

cache_mode = not args.prev_final_fasta        is None \
             and not args.prev_primers_outdir is None
if cache_mode:
    prev_fasta_fpath = os.path.abspath(args.prev_final_fasta)
    prev_primers_outdpath = os.path.abspath(args.prev_primers_outdir)
else:
    prev_fasta_fpath = None
    prev_primers_outdpath = None
# end if


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

# Primary check of "cached" data
if cache_mode:
    if not os.path.exists(prev_fasta_fpath):
        print('Error: file `{}` does not exist'.format(prev_fasta_fpath))
        sys.exit(1)
    # end def
    if not os.path.isdir(prev_primers_outdpath):
        print('Error: directory `{}` does not exist'.format(prev_primers_outdpath))
        sys.exit(1)
    # end def
# end def


print(fasta_fpath)
print(mfeprimer_fpath)
if cache_mode:
    print(prev_fasta_fpath)
    print(prev_primers_outdpath)
# end def
print()

# Configure k-mer for mfeprimer.
# Additionaly, it it minimal number of perfectly matching 3'-terminal primer bases.
K_MER_SIZE = 8


def prepare_pcr_template(seq_str):
    tmp_fasta = os.path.join(tmp_dirpath, 'tmpQ.fasta')

    # Write current sequence to fasta file. This will be input for mfeprimer.
    with open(tmp_fasta, 'wt') as tmp_fasta_file:
        tmp_fasta_file.write(f'>query_sequence\n{seq_str}\n')
    # end with

    # Index input sequence for mfeprimer
    index_fasta_for_mfeprimer(tmp_fasta)

    return tmp_fasta
# end def


def index_fasta_for_mfeprimer(fasta_fpath):
    index_cmd = f'{mfeprimer_fpath} index -i {fasta_fpath} -k {K_MER_SIZE} -c 2'
    os.system(index_cmd)
# end def


def simulate_pcr_for_single_template(template_fpath, primers_fpath):

    tmp_out_base = os.path.join(tmp_dirpath, 'tmpOUT')

    cmd = ' '.join(
        [
            mfeprimer_fpath, 'spec',
            f'--misEnd {K_MER_SIZE}',
            f'-k {K_MER_SIZE}',
            '-c 2',
            f'-i {primers_fpath}',
            f'-d {template_fpath}',
            # JSON mode remnant
            # '-j',
            # f'-o {tmp_out_base}'
        ]
    )
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, encoding='utf-8')
    stdout, stderr = pipe.communicate()
    if pipe.returncode != 0:
        print(f'MFEprimer exited with an error (exit code {pipe.returncode}):')
        print(stderr)
        sys.exit(1)
    # end if

    # JSON mode remnant
    # Clean
    # os.unlink(tmp_out_base)

    # JSON mode remnant
    # out_json_fpath = f'{tmp_out_base}.json'
    return stdout
# end def

# JSON mode remnant
# def parse_pcr_json_result(json_fpath):
#     # Read mfeprimer's output
#     mfe_json = json.loads(open(json_fpath, 'rt').read())

#     # Get list of possible amplicons (products)
#     amp_list = mfe_json['AmpList']

#     if not amp_list is None:
#         for amp in amp_list:
#             # Collect all appropriate data about primer annealing
#             product_size = amp['P']['Size']
#             ppc = amp['PPC'] # what's this?

#             f_size = amp['F']['Size'] # size of forw primer
#             f_start = amp['F']['Start'] # start pos of forw primer annealing
#             f_end = amp['F']['End'] # start pos of forw primer annealing
#             f_tm = amp['F']['Tm'] # melting temperature for forw primer
#             f_dg = amp['F']['Dg'] # free energy for forw primer
#             f_bind_len = len(amp['F']['Sseq']) # product sequence
#             f_pident = amp['F']['Aseq'].count(':') / len(amp['F']['Aseq']) # pident for forw primer
#             f_cover = f_bind_len / f_size # coverage for forw primer

#             r_size = amp['R']['Size'] # size of rev primer
#             r_start = amp['R']['Start'] # start pos of rev primer annealing
#             r_end = amp['R']['End'] # start pos of rev primer annealing
#             r_tm = amp['R']['Tm'] # melting temperature for rev primer
#             r_dg = amp['R']['Dg'] # free energy for rev primer
#             r_bind_len = len(amp['R']['Sseq']) # product sequence
#             r_pident = amp['R']['Aseq'].count(':') / len(amp['R']['Aseq']) # pident for rev primer
#             r_cover = r_bind_len / r_size # coverage for rev primer

#             output_values = list(
#                 map(
#                     str,
#                     [
#                         product_size, ppc,
#                         f_size, f_start, f_end, f_tm, f_dg, f_bind_len, f_pident, f_cover,
#                         r_size, r_start, r_end, r_tm, r_dg, r_bind_len, r_pident, r_cover,
#                     ]
#                 )
#             )

#             yield output_values
#         # end for
#     # end if

#     # Clean
#     os.unlink(json_fpath)

#     return
# # end def


def parse_pcr_plain_result(plain_text_str):
    num_ampl_reobj = re.search(
        r'Descriptions of \[ ([0-9]+) \] potential amplicons',
        plain_text_str
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

    for i_ampl in range(1, num_amplicons+1):
        pattern = PATTERN_DRAFT.format(i_ampl)
        re_obj = re.search(pattern, plain_text_str)

        if re_obj is None:
            print('ERROR')
            print('Cannot parse MFEPrimer stdout: cannot find PCR results')
            sys.exit(1)
        # end if

        output_values = [
            # TODO: remove
            # re_obj.group(1),  # t_start
            # re_obj.group(2),  # t_end
            re_obj.group(1),  # product_size
            re_obj.group(2),  # f_tm
            re_obj.group(3),  # f_dg
            re_obj.group(4),  # f_start
            re_obj.group(5),  # f_end
            re_obj.group(6),  # r_tm
            re_obj.group(7),  # r_dg
            re_obj.group(8), # r_start
            re_obj.group(9), # r_end
        ]

        yield output_values
    # end for
    return
# end def


def write_output_for_dedup_seq(output_row, seq_id_list, outfile):
    for seq_id in seq_id_list:
        # Get assembly ID of current seq_id
        asm_acc = parse_asm_acc(seq_id)
        # Add Assembly ID and SeqID to output_row
        output_row = list(map(str, output_row))
        complete_output_row = [asm_acc, seq_id] + output_row
        # Write output line
        outfile.write('{}\n'.format('\t'.join(complete_output_row)))
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


def read_cached_data(prev_fasta_fpath, prev_primers_outdpath):
    global PARTIAL_OUT_COLNAMES
    cached_seq_records = tuple(SeqIO.parse(prev_fasta_fpath, 'fasta'))
    uniq_cached_seqs = set(
        map(
            lambda r: str(r.seq),
            cached_seq_records
        )
    )
    print('Found {:,} unique cached sequences'.format(len(uniq_cached_seqs)))

    cache_dict = {
        seq: dict() for seq in uniq_cached_seqs
    }

    seq_2_cached_seqID = {
        seq: None for seq in uniq_cached_seqs
    }
    for r in cached_seq_records:
        seq = str(r.seq)
        if seq_2_cached_seqID[seq] is None:
            seq_2_cached_seqID[seq] = r.id
        # end if
    # end for
    seqIDs_for_cache_finding = set(seq_2_cached_seqID.values())

    del cached_seq_records

    print('Collecting cached data for primer pairs:')

    for i, (nameF, nameR) in enumerate(primer_pairs):
        primer_pair_key = make_primer_pair_key(nameF, nameR)

        print(
            '{} -- Start collecting for pair {}/{}: {}' \
                .format(get_time(), i+1, len(primer_pairs), primer_pair_key)
        )

        cached_fpath = os.path.join(
            prev_primers_outdpath,
            '{}.tsv'.format(primer_pair_key)
        )
        if not os.path.exists(cached_fpath):
            print('Error cached file does not exist: `{}`'.format(cached_fpath))
            sys.exit(1)
        # end if
        cached_df = pd.read_csv(cached_fpath, sep='\t')
        reduced_cached_df = cached_df.query('seqID in @seqIDs_for_cache_finding')
        del cached_df

        for cached_seq, cached_seqID in seq_2_cached_seqID.items():
            df_for_curr_seq = reduced_cached_df[
                reduced_cached_df['seqID'] == cached_seqID
            ][PARTIAL_OUT_COLNAMES]

            cache_dict[cached_seq][primer_pair_key] = [
                list(row) for _, row in df_for_curr_seq.iterrows()
            ]
        # end for
    # end for

    return cache_dict, uniq_cached_seqs
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
with open(primer_pairs_fpath, 'rt') as infile:
    primer_pairs = json.load(infile)
# end with


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

# Columns for output files
# JSON mode remnant
# OUT_COLNAMES = [
#     'asm_acc', 'seqID', 'product_size', 'ppc',
#     'f_size', 'f_start', 'f_end', 'f_tm', 'f_dg', 'f_bind_len', 'f_pident', 'f_cover',
#     'r_size', 'r_start', 'r_end', 'r_tm', 'r_dg', 'r_bind_len', 'r_pident', 'r_cover',
# ]
OUT_COLNAMES = [
    'asm_acc', 'seqID',
    # TODO: remove
    # 't_start', 't_end',
    'product_size', 
    'f_tm', 'f_dg', 'f_start', 'f_end', 
    'r_tm', 'r_dg', 'r_start', 'r_end', 
]

PARTIAL_OUT_COLNAMES = OUT_COLNAMES[2:]

# A regex pattern for searching PCR results parameters
PATTERN_DRAFT = r'''Amp {}: .+ \+ .+ ==> .+:[0-9]+-[0-9]+[ ]*
[ ]*
[ ]*Size = ([0-9]+) bp, GC content = [0-9\.]+%, Tm = [0-9\.-]+ .C, Ta = [0-9\.-]+ .C[ ]*
[ ]*F: Tm = ([0-9\.-]+) .C, Delta G = ([0-9\.-]+) kcal/mol, Start = ([0-9]+), End = ([0-9]+)[ ]*
[ ]*R: Tm = ([0-9\.-]+) .C, Delta G = ([0-9\.-]+) kcal/mol, Start = ([0-9]+), End = ([0-9]+)'''


print('{} -- Creating auxiliary data structures...'.format(get_time()))

# Read sequences
seq_records = tuple(SeqIO.parse(fasta_fpath, 'fasta'))
print('{} -- 1/4'.format(get_time()))

# Create output files for all primer pairs.
# And create temporary fasta files of primer pairs for mfeprimer.
primer_pair_fasta_dict = dict()

for nameF, nameR in primer_pairs:

    primer_pair_key = make_primer_pair_key(nameF, nameR)

    # Write header for output TSV files
    outfpath = primer_pair_key_2_outfpath(outdir_path, primer_pair_key)
    with open(outfpath, 'w') as outfile:
        outfile.write('{}\n'.format('\t'.join(OUT_COLNAMES)))
    # end with

    curr_primer_fpath = os.path.join(tmp_primers_dpath, f'{primer_pair_key}.fasta')
    # Write fasta data of current primer pair
    with open(curr_primer_fpath, 'w') as tmp_primers_file:
        tmp_primers_file.write(f'>{nameF}\n{primers[nameF]}\n>{nameR}\n{primers[nameR]}\n')
    # end with

    primer_pair_fasta_dict[primer_pair_key] = curr_primer_fpath
# end for
print('{} -- 2/4'.format(get_time()))


# De-replicate input sequences
uniq_seqs = set(
    map(lambda r: str(r.seq), seq_records)
)
uniq_seq_records = {
    seq: list() for seq in uniq_seqs
}
del uniq_seqs
for r in seq_records:
    uniq_seq_records[str(r.seq)].append(r.id)
# end for
print('{} -- 3/4'.format(get_time()))
del seq_records

if cache_mode:
    print('Reading cached data...')
    cache_dict, cached_seqs = read_cached_data(prev_fasta_fpath, prev_primers_outdpath)
else:
    cache_dict, cached_seqs = None, set()
# end if
print('{} -- 4/4'.format(get_time()))


# Count sequences for status bar
num_seqs = len(uniq_seq_records)
print('done\n')


# == Proceed ==
for i, (seq, curr_seq_id_list) in enumerate(uniq_seq_records.items()):
    print(
        '\r{} -- Doing deduplicated sequence {}/{}' \
            .format(get_time(), i+1, num_seqs),
        end=' '
    )

    clean_tmp_dir()

    if not seq in cached_seqs:
        template_fpath = prepare_pcr_template(seq)
    # end if

    # For each primer pair, check if they can anneal
    for nameF, nameR in primer_pairs:

        primer_pair_key = make_primer_pair_key(nameF, nameR)

        if cache_mode and seq in cached_seqs:
            output_rows_list = cache_dict[seq][primer_pair_key]
        else:
            # Get path to fasta file of current primer pair
            tmp_primers_fpath = primer_pair_fasta_dict[primer_pair_key]

            # Run mfeprimer
            tmp_out = simulate_pcr_for_single_template(
                template_fpath,
                tmp_primers_fpath
            )

            # JSON mode remnant
            # output_rows_list = parse_pcr_json_result(tmp_out)
            output_rows_list = parse_pcr_plain_result(tmp_out)
        # end if

        # Form path to current output file (corresponding to current primer pair)
        outfpath = primer_pair_key_2_outfpath(outdir_path, primer_pair_key)
        with open(outfpath, 'at') as outfile:
            for output_row in output_rows_list:
                # Write the result for every sequence identical to `seq`
                write_output_for_dedup_seq(output_row, curr_seq_id_list, outfile)
            # end def
        # end with
    # end for
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
