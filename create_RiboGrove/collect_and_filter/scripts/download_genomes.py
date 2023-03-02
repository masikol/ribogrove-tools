#!/usr/bin/env python3

# TODO: add description


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import time
import argparse

import pandas as pd

import src.rg_tools_IO as rgIO
from src.GenomeDownloader import GenomeDownloader, DownloadStatus
from src.rg_tools_time import print_time, get_time


# == Parse arguments ==

parser = argparse.ArgumentParser()

parser.add_argument(
    '-i',
    '--asm-sum',
    help="""TODO: add help""",
    required=True
)

parser.add_argument(
    '-o',
    '--outdir',
    help='output directory that will contain downloaded gbff.gz files',
    required=True
)

parser.add_argument(
    '-l',
    '--log-file',
    help='log file to track if all genomes are successfully downloaded',
    required=True
)


args = parser.parse_args()


asm_sum_fpath = os.path.realpath(args.asm_sum)
outdir = os.path.realpath(args.outdir)
log_fpath = os.path.realpath(args.log_file)


# Check existance of input file -c/--gi-2-acc-fpath
if not os.path.exists(asm_sum_fpath):
    print(f'Error: file `{asm_sum_fpath}` does not exist!')
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
del some_dir

print(asm_sum_fpath)
print()


def download_genomes(ass_sum_df, outdir, log_fpath):
    total_genome_count = ass_sum_df.shape[0]
    downloaded_count, already_here_count, failed_count = 0, 0, 0

    sys.stdout.write('{} -- 0/{:,} genomes done: 0 already here, 0 downloaded, 0 failed' \
        .format(get_time(), total_genome_count)
    )
    sys.stdout.flush()

    with open(log_fpath, 'wt') as logfile:
        for i, row in ass_sum_df.iterrows():
            asm_acc = row['asm_acc']

            downloader = GenomeDownloader(row, outdir)
            download_status = downloader.try_download()
            status_code = download_status.status_code

            if status_code == DownloadStatus.ALREADY_HERE:
                already_here_count += 1
                logfile.write(
                    '{} -- {}: genome is already here\n'.format(get_time(), asm_acc)
                )
            elif status_code == DownloadStatus.DOWNLOADED:
                downloaded_count += 1
                logfile.write(
                    '{} -- {}: ok, downloaded\n'.format(get_time(), asm_acc)
                )
            elif status_code == DownloadStatus.FAILED:
                failed_count += 1
                err_msg = download_status.err_msg.replace('\n', ' ')
                logfile.write(
                    '{} -- {}: failed. Error: {}\n' \
                        .format(get_time(), asm_acc, err_msg)
                )
            else:
                raise ValueError(
                    'Invalid GenomeDownloader.download status code: `{}`' \
                        .format(status_code)
                )
            # end if

            sys.stdout.write(
                '\r{} -- {:,}/{:,} genomes done: {:,} already here, {:,} downloaded, {:,} failed' \
                    .format(get_time(), i+1, total_genome_count, already_here_count, downloaded_count, failed_count)
            )
            sys.stdout.flush()
            time.sleep(1)
        # end for
    # end with

    sys.stdout.write('\n')
    sys.stdout.flush()
    return downloaded_count, already_here_count, failed_count
# end def


# == Proceed ==

ass_sum_df = rgIO.read_ass_sum_file(asm_sum_fpath)
downloaded_count, already_here_count, failed_count = download_genomes(
    ass_sum_df,
    outdir,
    log_fpath
)


total_genome_count = ass_sum_df.shape[0]

print('\n{} -- Completed!'.format(get_time()))
print('  {:,}/{:,} genomes were actually downloaded.'.format(downloaded_count, total_genome_count))
print('  {:,}/{:,} other genomes have been already downloaded earlier.'.format(already_here_count, total_genome_count))
print('  {:,}/{:,} genomes failed to get downloaded (see the log file for details).\n'.format(failed_count, total_genome_count))
print(outdir)
print(log_fpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')

if failed_count == 0:
    sys.exit(0)
else:
    sys.exit(1)
# end if
