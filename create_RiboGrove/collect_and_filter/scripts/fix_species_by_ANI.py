#!/usr/bin/env python3

# The script renames the original taxid column of an Assembly Summary table to raw_taxid
#   and adds a new column also named taxid, but these new taxids are ANI-corrected,
#   i.e. if the original (submitted) species identification is incorrect
#   and the genome has ANI>ANI_THRESHOLD and query_coverage>Q_COV_THRESHOLD when comparing to a type species i,
#   then the genome will acquire taxid of the species i.
# Thus many sp. genomes will acquire species names.


## Command line arguments

### Input files:
# 1. `-i / --in-asm-sum` -- an assembly summary file after the 4th step of filtering.
#   Mandatory.
# 2. `--in-raw-asm-sum` -- the raw, unfiltered assembly summary for taxid matching.
#   Mandatory.
# 3. `-g / --genomes-dir` -- a directory where the downlaoded genomes are located.
#   It is the output of the script `download_genomes.py`.
#   Mandatory.
# 3. `-m / --replicon-map` -- a replicon map file.
#   This is the output of the script `make_replicon_map.py`
#   Mandatory.

### Output files:

# 1. `-o / out-asm-sum` -- an amended assembly summary file.
#   Mandatory.


import os
import sys
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
    '-i',
    '--in-asm-sum',
    help='an assembly summary file after the 4th step of filtering',
    required=True
)

parser.add_argument(
    '--in-raw-asm-sum',
    help='the raw assembly summary file',
    required=True
)

parser.add_argument(
    '-g',
    '--genomes-dir',
    help='directory that contains downloaded gbff.gz files',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--out-asm-sum',
    help='a corrected assembly summary file',
    required=True
)

args = parser.parse_args()

infpath = os.path.realpath(args.in_asm_sum)
raw_asm_sum_fpath = os.path.realpath(args.in_raw_asm_sum)
genomes_dirpath = os.path.realpath(args.genomes_dir)
outfpath = os.path.realpath(args.out_asm_sum)

for f in (infpath, raw_asm_sum_fpath,):
    if not os.path.isfile(f):
        print('Error: file `{}` does not exist'.format(f))
        sys.exit(1)
    # end if
# end for
if not os.path.isdir(genomes_dirpath):
    print('Error: directory `{}` does not exist'.format(genomes_dirpath))
    sys.exit(1)
# end if


# == Import them now ==

import re
import json
import subprocess as sp
from pathlib import Path

import polars as pl

import src.rg_tools_IO as rgIO
from src.file_navigation import get_ani_report_fpath, \
                                get_ani_report_not_found_fpath


# README
# https://www.ncbi.nlm.nih.gov/datasets/docs/v2/data-processing/policies-annotation/quality/ani/
# https://pmc.ncbi.nlm.nih.gov/articles/PMC549018/
ANI_THRESHOLD = 95.0
Q_COV_THRESHOLD = 10.0

SUBM_ORGANISM_RE = re.compile(
    r'^Submitted organism:\s*(.*?)\s*\(taxid\s*=\s*(\d+)\s*,\s*rank\s*=\s*([^,)]+)',
    re.MULTILINE
)
BEST_MATCH_RE = re.compile(
    r'^Best match:\s*(.*?)\s*\(taxid\s*=\s*(\d+)\s*,\s*rank\s*=\s*([^,)]+)',
    re.MULTILINE
)
ANI_ROW_RE = re.compile(
    r'^\s*'
    r'(?P<ANI>\d+(?:\.\d+)?)\s+'
    r'\(\s*(?P<qcov>\d*(?:\.\d+)?)\s+\d*(?:\.\d+)?\s*\)\s+'
    r'\d+\s+\d+\s+(?:\d+\s+)?\d+\s+'
    r'(?:(?P<Flg>[CET])\s+)?'
    r'(?P<Organism>[A-Z].+) \('
    r'(?P<s_asm_acc>[GCAF_\.0-9]+),\s+'
    r'.+'
    r'\s*$',
    re.MULTILINE,
)


def make_ANI_renaming_dict(in_asm_sum_df, raw_asm_sum_df, genomes_dirpath) -> dict:
    ANI_renaming_dict = dict()

    for asm_acc in sorted(frozenset(in_asm_sum_df['asm_acc'])):
        ani_report_not_found_fpath = get_ani_report_not_found_fpath(
            asm_acc,
            genomes_dirpath
        )
        no_ani_report = os.path.isfile(ani_report_not_found_fpath)
        if no_ani_report:
            continue
        # end if

        ani_report_fpath = get_ani_report_fpath(asm_acc, genomes_dirpath)
        ani_df, subm_taxid, subm_name, new_taxid, new_tax_name = parse_ani_report(ani_report_fpath)

        ani_df = ani_df.filter(
            (pl.col('ANI') > ANI_THRESHOLD) & (pl.col('qcov') > Q_COV_THRESHOLD)
        ).sort(by='ANI', descending=True)
        if ani_df.height == 0:
            continue
        # end if

        new_tax_name_from_table = ani_df[0, 'Organism']
        if get_species(new_tax_name_from_table) != get_species(new_tax_name):
            # Fallback, if something weird happen, e.g.:
            # GCA_002945325.1
            # ANI     (Coverages)   NewSeq   CntmSeq  Assembly  Flg Organism  (assembly_accession, assembly_name)
            # ------- ------------- -------- -------- --------- --- --------------------------------------------------------------------
            #  81.396 ( 22.1  22.0)     9946     9946     17968     Methanococcus vannielii SB (GCA_000017165.1, ASM1716v1)
            #  95.633 ( 89.6  87.7)  1426889  1426889  26149618     Methanococcus maripaludis (GCA_017875815.1, ASM1787581v1)
            #  95.633 ( 89.0  88.0)   109785   109785  25036278     Methanococcus maripaludis (GCA_016907005.1, ASM1690700v1)
            #  82.541 (  2.8   2.7)        0        0  35392088     Methanococcus voltae PS (GCA_024807035.1, ASM2480703v1)
            #  80.217 (  1.8   1.9)        0        0    599208     Methanothermococcus thermolithotrophicus DSM 2095 (GCA_000376965.1, ASM37696v1)
            gb_asm_acc = ani_df[0, 's_asm_acc']
            new_taxid, new_tax_name = get_taxonomy_by_gb_asm_acc(gb_asm_acc, raw_asm_sum_df)
        # end if

        if new_taxid == subm_taxid:
            continue
        # end if
        if get_species(new_tax_name) == get_species(subm_name):
            continue
        # end if

        print('{} ({}, taxid {}) will be renamed to {}, taxid {}'.format(
            asm_acc,
            subm_name,
            subm_taxid,
            new_tax_name,
            new_taxid
        ))
        ANI_renaming_dict[asm_acc] = new_taxid
    # end for

    return ANI_renaming_dict
# end def


def parse_ani_report(path):
    text = Path(path).read_text()

    subm_taxid = None
    subm_name = None
    best_taxid = None
    best_name = None

    subm_organism = SUBM_ORGANISM_RE.search(text)
    if subm_organism is None:
        raise ValueError('Cannot find Submitter organism line in `{}`'.format(path))
    # end if
    subm_name = subm_organism.group(1).strip()
    subm_taxid = int(subm_organism.group(2))

    best_match = BEST_MATCH_RE.search(text)
    if best_match and best_match.group(3).strip().lower() == 'species':
        best_name = best_match.group(1).strip()
        best_taxid = int(best_match.group(2))
    # end if

    rows = [
        (
            float(match.group('ANI')),
            float(match.group('qcov')),
            match.group('Flg'),
            match.group('Organism'),
            match.group('s_asm_acc'),
        )
        for match in ANI_ROW_RE.finditer(text)
    ]

    schema = {
        'ANI': pl.Float64,
        'qcov': pl.Float64,
        'Flg': pl.Utf8,
        'Organism': pl.Utf8,
        's_asm_acc': pl.Utf8,
    }

    ani_df = pl.DataFrame(rows, schema=schema, orient='row')

    # Remove "contaminant rows"
    ani_df = (
        ani_df
        .with_columns(pl.col('Flg').replace({None: '-'}))
        .filter(pl.col('Flg') != 'C')
    )

    return ani_df, subm_taxid, subm_name, best_taxid, best_name
# end def


def get_taxonomy_by_gb_asm_acc(gb_asm_acc, raw_asm_sum_df):
    result_df = raw_asm_sum_df.filter(
        pl.col('gbrs_paired_asm') == gb_asm_acc
    )

    taxid, organism_name = None, None

    if result_df.height != 0:
        taxid = result_df[0,'taxid']
        organism_name = result_df[0,'organism_name']
    else:
        taxid, organism_name = get_taxonomy_via_ncbi_datasets(gb_asm_acc)
        print('# Got [{}, {}] for {} via ncbi datasets'.format(
            taxid,
            organism_name,
            gb_asm_acc
        ))
    # end if

    return taxid, organism_name
# end def

def get_taxonomy_via_ncbi_datasets(asm_acc):
    cmd = ' '.join([
        'datasets',
        'summary',
        'genome',
        'accession',
        asm_acc
    ])
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE, text=True, encoding='utf-8')
    out, err = pipe.communicate()
    if pipe.returncode != 0 or out.strip() == 'null':
        print(
            'Error: canot find assembly info for {} via ncbi datasets cli'.format(asm_acc),
            file=sys.stderr
        )
        print(err, file=sys.stderr)
        sys.exit(1)
    # end if

    try:
        data = json.loads(out)
    except json.decoder.JSONDecodeError as err:
        print(
            'Error: canot decode json requested for {}'.format(asm_acc),
            file=sys.stderr
        )
        print(err, file=sys.stderr)
        sys.exit(1)
    # end try

    try:
        org_info = data['reports'][0]['assembly_info']['biosample']['description']['organism']
        taxid = org_info['tax_id']
    except KeyError:
        return None, None
    # end try
    try:
        organism_name = org_info['organism_name']
    except KeyError:
        return taxid, 'NA'
    # end try

    return taxid, organism_name
# end def

def get_species(tax_name) -> str:
    return tax_name.split(' ')[:2]
# end def



# == Proceed ==

# Read input
in_asm_sum_df = rgIO.read_ass_sum_file(infpath)

raw_asm_sum_df = rgIO.read_ass_sum_file(raw_asm_sum_fpath, raw_summary=True)

# Make a dict that would determine how we will replace taxids
ANI_renaming_dict = make_ANI_renaming_dict(
    in_asm_sum_df,
    raw_asm_sum_df,
    genomes_dirpath
)

# Replace taxids
out_asm_sum_df = (
    in_asm_sum_df
    .clone()
    .rename({'taxid': 'raw_taxid'})
    .with_columns(
        pl.col('asm_acc').alias('taxid')
    )
    .with_columns(
        pl.col('taxid').replace_strict(
            ANI_renaming_dict,
            default=pl.col('raw_taxid')
        )
    )
)

out_asm_sum_df.write_csv(
    outfpath,
    separator='\t',
    null_value='NA',
    include_header=True,
    compression='gzip'
)


print('\n{} -- Completed!'.format(get_time()))
print(outfpath)
print(
    '\n|=== {} EXITTING SCRIPT `{}` ===|\n' \
    .format(
        get_time(), os.path.basename(__file__)
    )
)

sys.exit(0)
