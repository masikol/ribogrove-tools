#!/usr/bin/env python3

# The script takes RefSeq GI numbers from an input file (outputted by
#   the script `assembly2gi_numbers.py`) and finds corresponding RefSeq ACCESSION.VERSIONs and titles.
#   Output file is a TSV file of 3 columns (`gi_number`, `acc`, `title`).
# Requires Internet connection.
# For example, the script takes GI number 2075061612 (https://www.ncbi.nlm.nih.gov/nuccore/2075061612)
#   and finds corresponding ACCESSION.VERSION (NZ_CP063062.1)
#   and title ("Chlamydia suis strain 111 Ry chromosome, complete genome") for this record.

## Command line arguments
### Input files:
# 1. `-i / --gi-file` -- input TSV file mapping Assembly IDs to RefSeq GI numbers. Mandatory.

### Output files:
# 1. `-o / --outfile` -- output TSV file. Mandatory.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import re
import sys
import time
# import argparse
import requests

import pandas as pd

# from Bio import Entrez
# Entrez.email = 'maximdeynonih@gmail.com'


# parser = argparse.ArgumentParser()

# parser.add_argument(
#     '-i',
#     '--gi-file',
#     help='TSV file (with header) with Assembly IDs and GI numbers separated by tabs',
#     required=True
# )

# parser.add_argument(
#     '-o',
#     '--outfile',
#     help='file mapping RefSeq GI numbers to corresponding ACCESSION.VERSION\'s and titles',
#     required=True
# )

# args = parser.parse_args()

# For convenience
# gi_fpath = os.path.realpath(args.gi_file)
# outfpath = os.path.realpath(args.outfile)


# Archaea
# infpath = '/mnt/1.5_drive_0/RiboGrove/RiboGrove_workdirs/10.216/archaea/gene_seqs/archaea_all_collected_stats.tsv'
# outfpath = '/mnt/1.5_drive_0/RiboGrove/RiboGrove_workdirs/transition_to_11.217/archaea/gene_seqs/archaea_all_collected_stats.tsv'

# Bacteria
infpath = '/mnt/1.5_drive_0/RiboGrove/RiboGrove_workdirs/10.216/bacteria/gene_seqs/bacteria_all_collected_stats.tsv'
outfpath = '/mnt/1.5_drive_0/RiboGrove/RiboGrove_workdirs/transition_to_11.217/bacteria/gene_seqs/bacteria_all_collected_stats.tsv'


# Check existance of input file
if not os.path.exists(infpath):
    print(f'Error: file `{infpath}` does not exist!')
    sys.exit(1)
# end if

if not os.path.isdir(os.path.dirname(outfpath)):
    try:
        os.makedirs(os.path.dirname(outfpath))
    except OSError as err:
        print(f'Error: cannot create directory `{os.path.dirname(outfpath)}`')
        sys.exit(1)
    # end try
# end if

print(infpath)
print()


# Read input
in_df = pd.read_csv(
    infpath,
    sep='\t'
)

max_chunk_size = 50
n_done_ids = 0



# == Proceed ==

in_ass_ids = tuple(
    set(
        (ass_id for ass_id in in_df['ass_id'])
    )
)
out_ass_ids  = [None] * len(in_ass_ids)
out_asm_accs = [None] * len(in_ass_ids)

ass_id_pattern = re.compile(
    r'DocumentSummary uid=\"([0-9]+)\"'
)
asm_acc_pattern = re.compile(
    r'AssemblyAccession\>(GCF_[0-9\.]+)\</AssemblyAccession'
)

print('\r0/{}'.format(len(in_ass_ids)), end=' ')

for i in range(0, len(in_ass_ids), max_chunk_size):

    curr_ass_ids = tuple(
        map(
            str,
            tuple(in_ass_ids[i : i + max_chunk_size])
        )
    )
    chunk_size = len(curr_ass_ids)

    url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=assembly&id={}'.format(
        ','.join(curr_ass_ids)
    )

    success = False
    for _ in range(5):
        response = requests.get(url)
        response_cotent = response.content.decode('utf-8')

        curr_ass_ids = re.findall(ass_id_pattern, response_cotent)
        curr_asm_accs = re.findall(asm_acc_pattern, response_cotent)
        if len(curr_ass_ids) == chunk_size and len(curr_ass_ids) == chunk_size:
            success = True
            break
        # end if
    # end for

    if not success:
        print('Error! No success!')
        print(url)
        sys.exit(1)
    # end if

    for j, k in zip(range(i, i+chunk_size), range(chunk_size)):
        out_ass_ids[j] = curr_ass_ids[k]
        out_asm_accs[j] = curr_asm_accs[k]
    # end for

    # Wait a bit: we don't want NCBI to ban us :)
    time.sleep(0.5)

    n_done_ids += chunk_size
    print('\r{}/{}'.format(n_done_ids, len(in_ass_ids)), end=' ')
# end for

print('\r{}/{}'.format(n_done_ids, len(in_ass_ids)))

ass_acm_df = pd.DataFrame(
    {
        'ass_id': tuple(map(int, out_ass_ids)),
        'asm_acc': out_asm_accs,
    }
)

out_df = in_df.merge(
    ass_acm_df,
    on='ass_id',
    how='outer'
).copy()

out_df = out_df[
    [
        'ass_id',
        'asm_acc',
        'gi_number',
        'acc',
        'title',
        'seq_start_truncation',
        'improper_16S_annotation',
        'topology',
        'num_genes',
        'min_len',
        'max_len',
        'mean_len',
        'median_len',
    ]
]


out_df.to_csv(
    outfpath,
    sep='\t',
    index=False,
    header=True,
    na_rep='NA'
)


print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
