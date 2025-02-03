#!/usr/bin/env python3

import os
import sys
import json
import gzip

import pandas as pd


infpath = sys.argv[1]
# /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_workdirs/22.228/archaea/aberrations_and_heterogeneity/per_base_entropy.tsv.gz
outfpath  = sys.argv[2]
# /mnt/cager-beast/m.2.2/RiboGrove/RiboGrove_workdirs/22.228/archaea/aberrations_and_heterogeneity/per_base_entropy_new.json



with gzip.open(infpath, 'rt') as input_handle:
    per_base_df = pd.read_csv(
        input_handle,
        sep='\t',
        usecols=['asm_acc', 'entropy',]
    )
# end with


per_base_dict = {asm_acc: list() for asm_acc in frozenset(per_base_df['asm_acc'])}
num_genomes = len(per_base_dict)

n_rows = per_base_df.shape[0]
step = 10000
next_report = step

print()
sys.stdout.write('\r0/{} (0%)'.format(n_rows))
sys.stdout.flush()
for i, row in per_base_df.iterrows():
    per_base_dict[row['asm_acc']].append(row['entropy'])
    if i >= next_report:
        percent_done = round(
            (i+1) / n_rows * 100,
            2
        )
        sys.stdout.write('\r{}/{} ({}%)'.format(i+1, n_rows, percent_done))
        sys.stdout.flush()
        next_report += step
    # end if
# end for
percent_done = round(
    (i+1) / n_rows * 100,
    2
)
sys.stdout.write('\r{}/{} ({}%)'.format(i+1, n_rows, percent_done))
sys.stdout.flush()
print()

# print()
# for i, asm_acc in enumerate(per_base_dict.keys()):
#     sys.stdout.write('\r{} {}/{}{}'.format(asm_acc, i+1, num_genomes, ' '*10))
#     sys.stdout.flush()
#     wrk_df = per_base_df[per_base_df['asm_acc'] == asm_acc]
#     per_base_dict[asm_acc] = wrk_df['entropy'].to_list()
# # end for
# del per_base_df, wrk_df
# print()

with open(outfpath, 'wt') as output_handle:
    json.dump(per_base_dict, output_handle)
# end with

print('Completed! Have fun!')
sys.exit(0)

