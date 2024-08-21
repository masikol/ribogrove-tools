#! /usr/bin/env python3


import json
import gzip

import pandas as pd


def tsv_gz_to_df(infpath):
    with gzip.open(infpath, 'rt') as infile:
        df = pd.read_csv(infile, sep='\t')
    # end with
    return df
# end def


def df_to_dict(per_base_df):
    return {
        asm_acc: extract_per_base_entropy_list(asm_acc, per_base_df)
        for asm_acc in set(per_base_df['asm_acc'])
    }
# end def


def extract_per_base_entropy_list(asm_acc, per_base_df):
    return list(
        per_base_df[per_base_df['asm_acc'] == asm_acc]['entropy']
    )
# end def


def dict_to_json_gz_file(per_base_dict, outfpath):
    with gzip.open(outfpath, 'wt') as outfile:
        json.dump(per_base_dict, outfile)
    # end with
# end def

def summarize_json_gz(infpath):
    with gzip.open(infpath, 'rt') as infile:
        per_base_dict = json.load(infile)
    # end with
    for asm_acc, entropy_list in per_base_dict.items():
        entropy_list = tuple(map(float, entropy_list))
    # end for
# end def


# infpath = '/mnt/data/Max/RiboGrove/RiboGrove_workdirs/18.224/archaea/aberrations_and_heterogeneity/per_base_entropy.tsv.gz'
# outfpath = '/home/cager/Downloads/tmp/per_base_entropy.json.gz'

infpath = '/mnt/data/Max/RiboGrove/RiboGrove_workdirs/18.224/bacteria/aberrations_and_heterogeneity/per_base_entropy.tsv.gz'
outfpath = '/home/cager/Downloads/tmp/bacteria_per_base_entropy.json.gz'

df = tsv_gz_to_df(infpath)
per_base_dict = df_to_dict(df)
dict_to_json_gz_file(per_base_dict, outfpath)
