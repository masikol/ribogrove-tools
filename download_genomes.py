#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# Script downloads genomes in GenBank format, and saves them all to local directory.
# Output files will be all gziped.
# The script downloads records using efetch utility (https://www.ncbi.nlm.nih.gov/books/NBK25497/)

# Input file (-i/--assm-acc-file) is output of script merge_assID2acc_and_remove_WGS.py.

# Output files are:
# 1. `.gbk.gz` files in directory `outdir`
# 2. log file `log_fpath`, which will contain information about errors.
#    If record is successfully downloaded, the script writes "ok" to correcponding line of log file.


import os
import sys
import time
import http
import gzip
import argparse
from functools import partial

import pandas as pd
from Bio import Entrez


Entrez.email = 'maximdeynonih@gmail.com'


# == Parse arguments ==

parser = argparse.ArgumentParser()

parser.add_argument(
    '-i',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '-o',
    '--outdir',
    help='output directory that will contain downloaded gbk.gz files',
    required=True
)

parser.add_argument(
    '-l',
    '--log-file',
    help='log file to track if all genomes are successfully downloaded',
    required=True
)


args = parser.parse_args()


assm_acc_fpath = os.path.realpath(args.assm_acc_file)
outdir = os.path.realpath(args.outdir)
log_fpath = os.path.realpath(args.log_file)


# Check existance of input file -c/--gi-2-acc-fpath
if not os.path.exists(assm_acc_fpath):
    print(f'Error: file `{assm_acc_fpath}` does not exist!')
    sys.exit(1)
# end if

for some_dir in (outdir, os.path.dirname(log_fpath)):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end for


def write_to_file(text: str, log_fpath: str):
    # Function for writingto log file
    if text != '':
        with open(log_fpath, 'a') as logfile:
            logfile.write(f'{text}\n')
        # end with
    # end if
# end def write_to_file

write_to_log_file = partial(write_to_file, log_fpath=log_fpath)

# Empty log file
with open(log_fpath, 'w') as _:
    pass
# end with


# Read input
accs = tuple(
    pd.read_csv(
        assm_acc_fpath,
        sep='\t'
    )['acc']
)


n_accs = len(accs)
max_errors = 3 # We will terminate if 3 errors occur in a row
downloaded = 0 # amount of downloaded records


# == Proceed ==

for i, acc in enumerate(accs):

    print(f'\rDoing {i+1}/{n_accs}: {acc}', end=' '*10)

    # Configure output file path (e.g. `NZ_CP063178.1.gbk.gz`)
    outfpath = os.path.join(outdir, f'{acc}.gbk.gz')

    # Skip this record if it is already downloaded
    if os.path.exists(outfpath):
        write_to_log_file(f'{acc} - ok (already exists)')
        continue
    # end if

    downloaded += 1

    # Request record
    errors = 0
    while errors < max_errors:

        try:
            gbk_handle = Entrez.efetch(
                db='nucleotide',
                id=acc,
                # rettype='gb',
                rettype='gbwithparts',
                retmode='text',
                linkname='refseq'
            )

            # Write gbk content to output file
            with gzip.open(outfpath, 'wt') as outfile:
                outfile.write(gbk_handle.read())
            # end with

            gbk_handle.close()

        except IOError as err:
            err_msg = f'Error: cannot fetch GenBank file for `{acc}`: {err}'
            print(f'\n{err_msg}')
            write_to_log_file(err_msg)
        except (RuntimeError, http.client.IncompleteRead):
            errors += 1
        else:
            break
        # end try

        # Wait a bit: we don't want NCBI to ban us :)
        time.sleep(0.4)
    # end while

    if errors >= max_errors:
        err_msg = f'Error: 3 times RuntimeError occured for `{acc}`'
        print(f'\n{err_msg}')
        write_to_log_file(err_msg)
    else:
        write_to_log_file(f'{acc} - ok')
    # end if

    print(f'\r{i+1}/{n_accs} done', end=' '*10)
# end for

print('\n\nCompleted!')
print(f'Number of actually downloaded genomes = {downloaded}')
print(outdir)
print(log_fpath)
