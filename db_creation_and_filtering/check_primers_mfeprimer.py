#!/usr/bin/env python3

# The script takes fasta file with sequences and checks if hardcoded PCR primers
#   can produce some product with input sequences as templates.
#   The script uses MFEprimer (https://www.mfeprimer.com/) tool for PCR simulation.

## Command line arguments
### Input files:
# 1. `-f / --fasta-seqs-file` -- an input fasta file of template sequences.
#   This is the file, which is the output of th script `drop_repeats.py`. Mandatory.
# 2. `-c / --categories-file` -- the per-gene categories file.
#   This is the file, which is the output of th script `assign_genome_categories.py`. Mandatory.

### Output files:
# 1. `-o / --outdir` -- an output directory, where output TSV files
#   for each primer pair will be located. Mandatory.

### Dependencies:
# 1. --mfeprimer -- an [MFEprimer](https://www.mfeprimer.com/) executable.
#   Mandatory.

### "Cached" files:
# 1. `--prev-final-fasta` -- a fasta file of final RiboGrove sequences
#    of the previous RiboGrove release. Optional.
# 2. `--prev-primers-outdir` -- a directory where results of this script
#   are stored, but for the previous RiboGrove release. Optional.


import os
import sys
import glob
import json
import shutil
import argparse
import subprocess as sp

from Bio import SeqIO
import pandas as pd


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--fasta-seqs-file',
    help='fasta file of SSU gene sequences',
    required=True
)


parser.add_argument(
    '-c',
    '--categories-file',
    help='TSV file (with header) containing per-gene categories',
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

args = parser.parse_args()


# For convenience
fasta_fpath = os.path.abspath(args.fasta_seqs_file)
categories_fpath = os.path.abspath(args.categories_file)
mfeprimer_fpath = os.path.abspath(args.mfeprimer)
outdir_path = os.path.abspath(args.outdir)

cached_data_provided = not args.prev_final_fasta    is None \
                   and not args.prev_primers_outdir is None
if cached_data_provided:
    prev_fasta_fpath = os.path.abspath(args.prev_final_fasta)
    prev_primers_outdpath = os.path.abspath(args.prev_primers_outdir)
else:
    prev_fasta_fpath = None
    prev_primers_outdpath = None
# end if


# Check existance of all input files and dependencies
for fpath in (fasta_fpath, categories_fpath, mfeprimer_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directory if needed
if not os.path.isdir(outdir_path):
    try:
        os.makedirs(outdir_path)
    except OSError as err:
        print(f'Error: cannot create directory `{outdir_path}`')
        sys.exit(1)
    # end try
# end if

# Check if mfeprimer executable is actually executable
if not os.access(mfeprimer_fpath, os.X_OK):
    print(f'Error: file `{mfeprimer_fpath}` is not executable!')
    sys.exit(1)
# end if

# Primary check of "cached" data
if cached_data_provided:
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
print(categories_fpath)
print(mfeprimer_fpath)
if cached_data_provided:
    print(prev_fasta_fpath)
    print(prev_primers_outdpath)
# end def
print()

# Configure k-mer for mfeprimer.
# Additionaly, it it minimal number of perfectly matching 3'-terminal primer bases.
K_MER_SIZE = 8


def make_primer_pair_key(nameF, nameR):
    return f'{nameF}-{nameR}'
# end def


def primer_pair_key_2_outfpath(outdir_path, primer_pair_key):
    return os.path.join(outdir_path, f'{primer_pair_key}.tsv')
# end def


def prepare_pcr_template(seq_str):
    tmp_fasta = os.path.join(tmp_dir, 'tmpQ.fasta')

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


    tmp_out_base = os.path.join(tmp_dir, 'tmpOUT')

    cmd = ' '.join(
        [
            mfeprimer_fpath, 'spec',
            f'--misEnd {K_MER_SIZE}',
            f'-k {K_MER_SIZE}',
            '-c 2',
            f'-i {primers_fpath}',
            f'-d {template_fpath}',
            '-j',
            f'-o {tmp_out_base}'
        ]
    )
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()
    if pipe.returncode != 0:
        print(f'MFEprimer exited with error (exit code {pipe.returncode}):')
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(1)
    # end if

    # Clean
    os.unlink(tmp_out_base)

    out_json_fpath = f'{tmp_out_base}.json'
    return out_json_fpath
# end def


def parse_pcr_json_result(json_fpath):
    # Read mfeprimer's output
    mfe_json = json.loads(open(json_fpath, 'rt').read())

    # Get list of possible amplicons (products)
    amp_list = mfe_json['AmpList']

    if not amp_list is None:
        for amp in amp_list:
            # Collect all appropriate data about primer annealing
            product_size = amp['P']['Size']
            ppc = amp['PPC'] # what's this?

            f_size = amp['F']['Size'] # size of forw primer
            f_start = amp['F']['Start'] # start pos of forw primer annealing
            f_end = amp['F']['End'] # start pos of forw primer annealing
            f_tm = amp['F']['Tm'] # melting temperature for forw primer
            f_dg = amp['F']['Dg'] # free energy for forw primer
            f_bind_len = len(amp['F']['Sseq']) # product sequence
            f_pident = amp['F']['Aseq'].count(':') / len(amp['F']['Aseq']) # pident for forw primer
            f_cover = f_bind_len / f_size # coverage for forw primer

            r_size = amp['R']['Size'] # size of rev primer
            r_start = amp['R']['Start'] # start pos of rev primer annealing
            r_end = amp['R']['End'] # start pos of rev primer annealing
            r_tm = amp['R']['Tm'] # melting temperature for rev primer
            r_dg = amp['R']['Dg'] # free energy for rev primer
            r_bind_len = len(amp['R']['Sseq']) # product sequence
            r_pident = amp['R']['Aseq'].count(':') / len(amp['R']['Aseq']) # pident for rev primer
            r_cover = r_bind_len / r_size # coverage for rev primer

            output_values = list(
                map(
                    str,
                    [
                        product_size, ppc,
                        f_size, f_start, f_end, f_tm, f_dg, f_bind_len, f_pident, f_cover,
                        r_size, r_start, r_end, r_tm, r_dg, r_bind_len, r_pident, r_cover,
                    ]
                )
            )

            yield output_values
        # end for
    # end if

    # Clean
    os.unlink(json_fpath)

    return
# end def


def write_output_for_dedup_seq(output_values, seq_id_list, outfpath):
    with open(outfpath, 'at') as outfile:
        for seq_id in curr_seq_id_list:
            # Get assembly ID of current seq_id
            ass_id = seqID_2_assID_dict[seq_id]
            # Add Assembly ID and SeqID to output_values
            output_values = list(map(str, output_values))
            complete_output_values = [str(ass_id), str(seq_id)] + output_values
            # Write output line
            outfile.write('{}\n'.format('\t'.join(complete_output_values)))
        # end for
    # end with
# end def


def clean_tmp_dir():
    # Remove all files from temporary directory
    for f in glob.iglob(f'{tmp_dir}/*'):
        if not os.path.isdir(f):
            os.unlink(f)
        # end if
    # end for
# end def


def read_cached_data(prev_fasta_fpath, prev_primers_outdpath):
    cached_seq_records = tuple(SeqIO.parse(prev_fasta_fpath, 'fasta'))
    uniq_cached_seqs = set(
        map(
            lambda r: str(r.seq),
            cached_seq_records
        )
    )

    cache_dict = {
        seq: dict() for seq in uniq_cached_seqs
    }

    seq_2_any_seqID = {
        seq: None for seq in uniq_cached_seqs
    }
    for r in cached_seq_records:
        seq = str(r.seq)
        if seq_2_any_seqID[seq] is None:
            seq_2_any_seqID[seq] = r.id
        # end if
    # end for

    del cached_seq_records

    print('Primer pairs:')

    for i, (nameF, nameR) in enumerate(primer_pairs):
        primer_pair_key = make_primer_pair_key(nameF, nameR)

        print('{}/{}: {}'.format(i+1, len(primer_pairs), primer_pair_key))

        cached_fpath = os.path.join(
            prev_primers_outdpath,
            '{}.tsv'.format(primer_pair_key)
        )
        if not os.path.exists(cached_fpath):
            print('Error cached file does not exist: `{}`'.format(cached_fpath))
            sys.exit(1)
        # end if
        cached_df = pd.read_csv(cached_fpath, sep='\t')

        necessary_seqIDs = set(seq_2_any_seqID.values())
        reduced_cached_df = cached_df.query('seqID in @necessary_seqIDs')
        del necessary_seqIDs, cached_df

        for cached_seq, any_seqID in seq_2_any_seqID.items():
            df_for_curr_seq = reduced_cached_df[
                reduced_cached_df['seqID'] == any_seqID
            ][partial_out_colnames]

            cache_dict[cached_seq][primer_pair_key] = list()
            for i, row in df_for_curr_seq.iterrows():
                values = list(row)
                cache_dict[cached_seq][primer_pair_key].append(values)
            # end for
        # end for
    # end for

    return cache_dict, uniq_cached_seqs
# end def


primers = {
    '27F': 'AGAGTTTGATYMTGGCTCAG',
    '338R': 'GCTGCCTCCCGTAGGAGT',
    '534R': 'ATTACCGCGGCTGCTGG',
    '341F': 'CCTACGGGNGGCWGCAG',
    '785R': 'GACTACHVGGGTATCTAATCC',
    '515F': 'GTGCCAGCMGCCGCGGTAA',
    '784F': 'AGGATTAGATACCCTGGTA',
    '1100R': 'AGGGTTGCGCTCGTTG',
    '806R': 'GGACTACHVGGGTWTCTAAT',
    '944R': 'GAATTAAACCACATGCTC',
    '939F': 'GAATTGACGGGGGCCCGCACAAG',
    '1115F': 'CAACGAGCGCAACCCT',
    '1193R': 'ACGTCATCCCCACCTTCC',
    '1378R': 'CGGTGTGTACAAGGCCCGGGAACG',
    '1492R': 'TACGGYTACCTTGTTACGACTT',
}

primer_pairs = [
    # Full
    ['27F', '1492R'],
    # Two
    # V1-V2
    ['27F', '338R'],
    # V3-V4
    ['341F', '785R'],
    # V4-V5
    ['515F', '944R'],
    # V5-V6
    ['784F', '1100R'],
    # V6-V7
    ['939F', '1193R'],
    # Three
    # V1-V3
    ['27F', '534R'],
    # V3-V5
    ['341F', '944R'],
    # V4-V6
    ['515F', '1100R'],
    # V5-V7
    ['784F', '1193R'],
    # V6-V8
    ['939F', '1378R'],
    # V7-V9
    ['1115F', '1492R'],
]


# Configure paths to temporary files
tmp_dir = os.path.join(outdir_path, 'tmp')
tmp_primers_dpath = os.path.join(tmp_dir, 'tmp_primers')


for some_dir in (tmp_dir, tmp_primers_dpath):
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
out_colnames = [
    'ass_id', 'seqID', 'product_size', 'ppc',
    'f_size', 'f_start', 'f_end', 'f_tm', 'f_dg', 'f_bind_len', 'f_pident', 'f_cover',
    'r_size', 'r_start', 'r_end', 'r_tm', 'r_dg', 'r_bind_len', 'r_pident', 'r_cover',
]

partial_out_colnames = out_colnames[2:]

# Create dictionary mapping seqIDs to Assembly IDs
print('Creating auxiliary data structures...')
seqID_2_assID_dict = {
    row['seqID']: row['ass_id'] \
    for _, row in pd.read_csv(categories_fpath, sep='\t').iterrows()
}
print('.')


# Create output files for all primer pairs.
# And create temporary fasta files of primer pairs for mfeprimer.
primer_pair_fasta_dict = dict()

for nameF, nameR in primer_pairs:

    primer_pair_key = make_primer_pair_key(nameF, nameR)

    # Write header for output TSV files
    outfpath = primer_pair_key_2_outfpath(outdir_path, primer_pair_key)
    with open(outfpath, 'w') as outfile:
        outfile.write('{}\n'.format('\t'.join(out_colnames)))
    # end with

    curr_primer_fpath = os.path.join(tmp_primers_dpath, f'{primer_pair_key}.fasta')
    # Write fasta data of current primer pair
    with open(curr_primer_fpath, 'w') as tmp_primers_file:
        tmp_primers_file.write(f'>{nameF}\n{primers[nameF]}\n>{nameR}\n{primers[nameR]}\n')
    # end with

    primer_pair_fasta_dict[primer_pair_key] = curr_primer_fpath
# end for
print('.')

# Read sequences
seq_records = tuple(SeqIO.parse(fasta_fpath, 'fasta'))
print('.')

# De-replicate input sequences
uniq_seqs = set(
    map(lambda r: str(r.seq), seq_records)
)
uniq_seq_records = {
    seq: list() for seq in uniq_seqs
}
for r in seq_records:
    seq = str(r.seq)
    uniq_seq_records[seq].append(r.id)
# end for
print('.')

del uniq_seqs
del seq_records

if cached_data_provided:
    print('Reading cached data...')
    cache_dict, cached_seqs = read_cached_data(prev_fasta_fpath, prev_primers_outdpath)
else:
    cache_dict, cached_seqs = None, set()
# end if
print('.')


# Count sequences for status bar
num_seqs = len(uniq_seq_records)
print('done\n')


# == Proceed ==
for i, (seq, curr_seq_id_list) in enumerate(uniq_seq_records.items()):
    print(f'\rDoing deduplicated sequence {i+1}/{num_seqs}', end=' ')

    clean_tmp_dir()

    if not seq in cached_seqs:
        template_fpath = prepare_pcr_template(seq)
    # end if

    # For each primer pair, check if they can anneal
    for nameF, nameR in primer_pairs:

        primer_pair_key = make_primer_pair_key(nameF, nameR)
        # Form path to current output file (corresponding to current primer pair)
        outfpath = primer_pair_key_2_outfpath(outdir_path, primer_pair_key)

        if seq in cached_seqs:
            output_values_collection = cache_dict[seq][primer_pair_key]
        else:
            # Get path to fasta file of current primer pair
            tmp_primers_fpath = primer_pair_fasta_dict[primer_pair_key]

            # Run mfeprimer
            tmp_out_json = simulate_pcr_for_single_template(
                template_fpath,
                tmp_primers_fpath
            )

            output_values_collection = parse_pcr_json_result(tmp_out_json)
        # end if

        for output_values in output_values_collection:
            # Write the result for every sequence identical to `seq`
            write_output_for_dedup_seq(output_values, curr_seq_id_list, outfpath)
        # end def
    # end for
# end for


# Remove temporary directory
try:
    shutil.rmtree(tmp_dir)
except OSError as err:
    print(f'Warning: Cannot delete temporary directory: `{tmp_dir}`')
    print(err)
    print()
# end try


print('\nCompleted!')
print(outdir_path)
