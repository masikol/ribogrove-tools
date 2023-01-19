#!/usr/bin/env python3

# This read_and_filter_fasta function reads a fasta file,
#   removes sequences having ids listed in one of the filter files,
#   removes sequences listed in the blacklist
#   but keeps sequences listed in the whitelist.

import gzip

from Bio import SeqIO


def read_and_filter_fasta(in_fasta_fpath,
                          filter_fpaths=[],
                          blacklist=set(),
                          whitelist=set()):

    seqIDs_to_rm = _get_seqIDs_to_rm(filter_fpaths, blacklist, whitelist)
    is_passing = lambda record: not record.id in seqIDs_to_rm

    if in_fasta_fpath.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    # end if

    with open_func(in_fasta_fpath, 'rt') as infile:
        seq_records = list(
            filter(
                is_passing,
                SeqIO.parse(infile, 'fasta')
            )
        )
    # end with

    return seq_records
# end def


def _get_seqIDs_to_rm(filter_fpaths, blacklist, whitelist):
    seqIDs_to_rm = set()
    for filter_fpath in filter_fpaths:
        with open(filter_fpath, 'rt') as infile:
            seqIDs_to_rm = seqIDs_to_rm | set(
                map(str.strip, infile.readlines())
            )
        # end with
    # end for

    seqIDs_to_rm = (seqIDs_to_rm | blacklist) - whitelist

    return seqIDs_to_rm
# end def
