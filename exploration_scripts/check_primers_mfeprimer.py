#!/usr/bin/env python3

# Script takes fasta file with sequences and checks if hardcoded primers
#   can produce some product with input sequences as templates.
# Script processes only sequences of category 1.
#
# Input files:
# 1. Input fasta file of template sequences (-f/--fasta-seqs-file)
# 2. Per-gene categories file (-c/--categories-file)
#
# Output files:
# TSV files for each primer pair in directory -o/--outdir
#
# Dependencies:
# 1. MFEprimer executable (https://www.mfeprimer.com/)

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
outdpath = os.path.abspath(args.outdir)


# Check existance of all input files and dependencies
for fpath in (fasta_fpath, categories_fpath, mfeprimer_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Create output directory if needed
if not os.path.isdir(outdpath):
    try:
        os.makedirs(outdpath)
    except OSError as err:
        print(f'Error: cannot create directory `{outdpath}`')
        sys.exit(1)
    # end try
# end if

# Check if mfeprimer executable is actually executable
if not os.access(mfeprimer_fpath, os.X_OK):
    print(f'Error: file `{mfeprimer_fpath}` is not executable!')
    sys.exit(1)
# end if


print(fasta_fpath)
print(categories_fpath)
print(mfeprimer_fpath)
print()


# Bacteria

primers = {
    # '27F': 'AGAGTTTGATCMTGGCTCAG',
    # '16S_BV2f': 'AGTGGCGGACGGGTGAGTAA',
    # 'CCF': 'CCAGACTCCTACGGGAGGCAGC',
    # '16SAmpliconPCRForwardPrimer': 'CCTACGGGNGGCWGCAG',
    # '336F-seq': 'ACGGYCCAGACTCCTACG',
    # '338R': 'GCTGCCTCCCGTAGGAGT',
    # '515F': 'GTGCCAGCMGCCGCGGTAA',
    # '522R-seq': 'TATTACCGCGGCTGCTGGCAC',
    # '16S_BV3r': 'CCGCGGCTGCTGGCAC',
    # '519R': 'GWATTACCGCGGCKGCTG',
    # '806R': 'GGACTACHVGGGTWTCTAAT',
    # '784F': 'AGGATTAGATACCCTGGTA',
    # '16SAmpliconPCRReversePrimer': 'GACTACHVGGGTATCTAATCC',
    # '907R': 'CCGYCAATTCMTTTRAGTTT',
    # '918F-seq': 'ACTCAAAKGAATTGACGGG',
    # '926R-seq': 'CCGTCAATTCATTTGAGTTT',
    # '968F': 'AACGCGAAGAACCTTAC',
    # 'CDR': 'CTTGTGCGGGCCCCCGTCAATTC',
    # '1061R': 'CRRCACGAGCTGACGAC',
    # '1100R': 'AGGGTTGCGCTCGTTG',
    # '1114F': 'GCAACGAGCGCAACCC',
    # '1193R': 'ACGTCATCCCCACCTTCC',
    # # '1237F': 'GGGCTACACACGYGCWAC',
    # '1391R': 'GACGGGCGGTGTGTRCA',
    # '1492Rl': 'GGTTACCTTGTTACGACTT',
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

# primer_pairs = [
#     # Full
#     ['27F', '1492Rl'],
#     # Two
#     # V1-V2
#     ['27F', '338R'],
#     # V2-V3
#     ['16S_BV2f', '16S_BV3r'],
#     ['16S_BV2f', '519R'],
#     ['16S_BV2f', '522R-seq'],
#     # V3-V4
#     ['16SAmpliconPCRForwardPrimer', '16SAmpliconPCRReversePrimer'],
#     ['16SAmpliconPCRForwardPrimer', '806R'],
#     ['CCF', '16SAmpliconPCRReversePrimer'],
#     ['CCF', '806R'],
#     ['336F-seq', '16SAmpliconPCRReversePrimer'],
#     ['336F-seq', '806R'],
#     # V4-V5
#     ['515F', '907R'],
#     ['515F', 'CDR'],
#     ['515F', '926R-seq'],
#     # V5-V6
#     ['784F', '1061R'],
#     ['784F', '1100R'],
#     # V6-V7
#     ['968F', '1193R'],
#     ['918F-seq', '1193R'],
#     # V7-V8
#     ['1114F', '1391R'],
#     # Three
#     # V1-V3
#     ['27F', '16S_BV3r'],
#     ['27F', '519R'],
#     ['27F', '522R-seq'],
#     # V2-V4
#     ['16S_BV2f', '16SAmpliconPCRReversePrimer'],
#     ['16S_BV2f', '806R'],
#     # V3-V5
#     ['16SAmpliconPCRForwardPrimer', '907R'],
#     ['16SAmpliconPCRForwardPrimer', 'CDR'],
#     ['CCF', '907R'],
#     ['CCF', 'CDR'],
#     ['336F-seq', '907R'],
#     ['336F-seq', 'CDR'],
#     ['16SAmpliconPCRForwardPrimer', '926R-seq'],
#     ['CCF', '926R-seq'],
#     # V4-V6
#     ['515F', '1061R'],
#     ['515F', '1100R'],
#     # V5-V7
#     ['784F', '1193R'],
#     # V6-V8
#     ['968F', '1391R'],
#     ['918F-seq', '1391R'],
#     # V7-V9
#     ['1114F', '1492Rl'],
# ]


# Configure paths to temporary files
tmp_dir = os.path.join(outdpath, 'tmp')
tmp_fasta = os.path.join(tmp_dir, 'tmpQ.fasta')
tmp_out_base = os.path.join(tmp_dir, 'tmpOUT')
tmp_out_json = f'{tmp_out_base}.json'
tmp_primers_dpath = os.path.join(tmp_dir, 'tmp_primers')

# Configure k-mer for mfeprimer.
# Additionaly, it it minimal number of perfectly matching 3'-terminal primer bases.
k_mer_size = 8

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
out_colnames = (
    'ass_id', 'seqID', 'product_size', 'ppc',
    'f_size', 'f_start', 'f_end', 'f_tm', 'f_dg', 'f_bind_len', 'f_pident', 'f_cover',
    'r_size', 'r_start', 'r_end', 'r_tm', 'r_dg', 'r_bind_len', 'r_pident', 'r_cover',
)

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

    # Kind of ID for primer pair
    primer_pair_key = f'{nameF}-{nameR}'

    # Write header for output TSV files
    outfpath = os.path.join(outdpath, f'{primer_pair_key}.tsv')
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

# Count sequences for status bar
num_seqs = len(seq_records)
print('done\n')


# == Proceed ==

for i, seq_record in enumerate(seq_records):

    print(f'\rDoing sequence {i+1}/{num_seqs}: {seq_record.id}', end=' ')

    # Get assembly ID of current sequence
    ass_id = seqID_2_assID_dict[seq_record.id]

    # Remove all files from temporary directory
    for f in glob.iglob(f'{tmp_dir}/*'):
        if not os.path.isdir(f):
            os.unlink(f)
        # end if
    # end for

    # Write current sequence to fasta file. This will be input for mfeprimer.
    with open(tmp_fasta, 'wt') as tmp_fasta_file:
        tmp_fasta_file.write(f'>{seq_record.id}\n{seq_record.seq}\n')
    # end with

    # Index input sequence for mfeprimer
    index_cmd = f'{mfeprimer_fpath} index -i {tmp_fasta} -k {k_mer_size} -c 2'
    os.system(index_cmd)

    # For feach primer pair, check if they can anneal
    for nameF, nameR in primer_pairs:

        # Kind of ID for primer pair
        primer_pair_key = f'{nameF}-{nameR}'

        # Form path to current output file (corresponding to current primer pair)
        outfpath = os.path.join(outdpath, f'{primer_pair_key}.tsv')

        # Get path to fasta file of current primer pair
        tmp_primers_fpath = primer_pair_fasta_dict[primer_pair_key]

        # Run mfeprimer
        cmd = f'{mfeprimer_fpath} spec --misEnd {k_mer_size} -k {k_mer_size} -c 2 -i {tmp_primers_fpath} -d {tmp_fasta} -j -o {tmp_out_base}'
        pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
        stdout_stderr = pipe.communicate()
        if pipe.returncode != 0:
            print(f'MFEprimer exited with error (exit code {pipe.returncode}):')
            print(stdout_stderr[1].decode('utf-8'))
            sys.exit(1)
        # end if

        # Read mfeprimer's output
        mfe_json = json.loads(open(tmp_out_json, 'rt').read())

        # Get list of possible amplicons (products)
        amp_list = mfe_json['AmpList']

        if not amp_list is None:
            with open(outfpath, 'at') as outfile:

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

                    # Form output values
                    out_values = map(
                        str,
                        (
                            ass_id, seq_record.id, product_size, ppc,
                            f_size, f_start, f_end, f_tm, f_dg, f_bind_len, f_pident, f_cover,
                            r_size, r_start, r_end, r_tm, r_dg, r_bind_len, r_pident, r_cover,
                        )
                    )

                    # Wrtire output line
                    outfile.write('{}\n'.format('\t'.join(out_values)))
                # end for
            # end with
        # end if

        # Remove temporary files for next run of mfeprimer
        os.unlink(tmp_out_base)
        os.unlink(tmp_out_json)
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
print(outdpath)
