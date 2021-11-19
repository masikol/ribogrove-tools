#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

# The script takes IDs of the NCBI Assembly database from an input file and finds
#   corresponding RefSeq GI numbers using elink utility (https://www.ncbi.nlm.nih.gov/books/NBK25497/).
#   The script writes result to the output TSV file of 2 columns (`ass_id`, `gi_numbers`),
#   thus mapping Assembly IDs to RefSeq GI numbers.
# Requires Internet connection.
# For example, the script takes Assembly ID 10601591 (https://www.ncbi.nlm.nih.gov/assembly/10601591/)
#   and finds corresponding RefSeq record [2075061612](https://www.ncbi.nlm.nih.gov/nuccore/2075061612).

## Command line arguments

### Input files:
# 1. `-i / --accs-file` -- input file of Assembly IDs (one per line). Mandatory.

###  Output files:
# 1. `-o / --outfile` -- output TSV file mapping Assembly IDs to RefSeq GI numbers. Mandatory.

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import time
import argparse
import http.client


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
    '--outfile',
    help='file mapping Assembly IDs to RefSeq GI numbers',
    required=True
)

args = parser.parse_args()


# For convenience
assm_accs_fpath = os.path.realpath(args.assm_acc_file)
outfpath = os.path.realpath(args.outfile)


# Check existance of input file
if not os.path.exists(assm_accs_fpath):
    print(f'Error: file `{assm_accs_fpath}` does not exist!')
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



def esummary_assembly(ass_id):

    # We will terminate if 3 errors occur in a row
    n_errors = 0

    # Request RefSeq GI number
    while n_errors < 3:
        try:
            handle = Entrez.esummary(
                db='assembly',
                id=ass_id
            )
            record = Entrez.read(handle)
            handle.close()
            break
        except (IOError, RuntimeError, http.client.HTTPException) as err:
            print('\n' + str(err))
            n_errors += 1
        # end try

        if len(record[0]['ERROR']) != 0:
            print('\n' + '\n'.join(record[0]['ERROR']))
            n_errors += 1
            print('Sleeping 30 seconds')
            time.sleep(30)
        # end if
    # end while

    if n_errors >= 3:
        print(f'\n3 errors (Assembly ID {ass_id})')
        print('Sleeping 60 seconds')
        time.sleep(60)
        return None
    # end if

    return record
# end def esummary_assembly


# Read all assembly IDs
ass_ids = set(
    map(
        lambda x: x.split('\t')[0],
        map(
            str.strip,
            open(assm_accs_fpath).readlines()[1:]
        )
    )
)

ass_ids = {1270421, 10960041}


# == Proceed ==

with open(outfpath, 'wt') as outfile:

    # Write header
    outfile.write('ass_id\tstatus\n')

    # Iterate over assembly IDs
    for i, ass_id in enumerate(ass_ids):

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' ')

        record = esummary_assembly(ass_id)
        if record is None:
            continue
        # end if

        ass_status = record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyStatus']
        # updated_uid = updated_record['DocumentSummarySet']['DocumentSummary'][0].attributes['uid']
        # init_ass_accession = record['DocumentSummarySet']['DocumentSummary'][0]['AssemblyAccession']

        outfile.write(f'{ass_id}\t{ass_status}\n')

        # Wait a bit: we don't want NCBI to ban us :)
        time.sleep(0.5)
    # end for
# end with

print('\nCompleted!')
print(outfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
