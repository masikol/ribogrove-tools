#!/usr/bin/env python3

# This file contains the function create_matrix, which calculates a frequency matrix
#   (a pandas DataFrame) for sequence logo creation using logomaker
#   https://logomaker.readthedocs.io/en/latest/.

# The function usage:
#
# matrix, coverage_array = create_matrix(fasta_fpath)


import os
import sys
import math
import operator
from functools import reduce

import numpy as np
import pandas as pd
from Bio import SeqIO


def _check_lengths(seq_records) -> int:
    # The function checks if all input sequences are of equal length (after MSA).

    distinct_lengths = set(
        map(
            len,
            seq_records
        )
    )

    if len(distinct_lengths) > 1:
        ValueError(f'Input fasta contains sequences of different lengths. Lengths: {distinct_lengths}')
    elif len(distinct_lengths) == 0:
        ValueError('Seems, input file does not contain fasta records.')
    # end if

    return next(iter(distinct_lengths))
# end def _check_lengths


def _extract_column(seqs, i):
    # The function extracts an alignment (MSA) column from MSA.

    return tuple(
        map(
            lambda x: x[i],
            seqs
        )
    )
# end def _extract_column


# def _calc_information(col_bases, base_vocabulary):

#     aln_column = list(
#         ''.join(col_bases).replace('-', '')
#     )

#     n_seqs = len(aln_column)

#     max_information = math.log(len(base_vocabulary), 2)

#     freq_dict = {base: aln_column.count(base) / n_seqs for base in base_vocabulary}

#     # Calculate information
#     # abs instead of minus in order not to allow "-0.0" values
#     information = max_information - abs(
#         reduce(
#             operator.add,
#             (
#                 freq * math.log(freq, 2) \
#                     for base, freq in filter(
#                         lambda x: x[0] in aln_column,
#                         freq_dict.items()
#                     )
#             )
#         )
#     )

#     information_dict = {base: 0.0 for base in base_vocabulary}

#     for base in information_dict.keys():
#         base_freq = freq_dict[base]
#         if base_freq > 1e-9:
#             information_dict[base] = information - ((-1) * base_freq * math.log(base_freq, 2))
#         # end if
#     # end for

#     coverage = n_seqs

#     return information_dict, coverage
# # end def _calc_information


def _calc_frequencies(col_bases, base_vocabulary):
    # The function calculates frequencies of bases.

    aln_column = list(
        ''.join(col_bases).replace('-', '')
    )

    n_seqs = len(aln_column)

    freq_dict = {base: aln_column.count(base) / n_seqs for base in base_vocabulary}

    coverage = n_seqs

    return freq_dict, coverage
# end def _calc_frequencies


def create_matrix(fasta_fpath):
    """
    The function creates a frequency matrix (a pandas DataFrame) for logomaker.
    The function usage:
        matrix, coverage_array = create_matrix(fasta_fpath)
    Where `fasta_fpath` is path of fasta file of sequences to make logo.
    """

    base_vocabulary = ['A', 'T', 'G', 'C']

    seq_records = tuple(SeqIO.parse(fasta_fpath, 'fasta'))
    seqs = tuple(
        map(
            lambda x: str(x.seq),
            seq_records
        )
    )
    logo_length = _check_lengths(seq_records)

    matrix_columns = {base: np.repeat(np.nan, logo_length) for base in base_vocabulary}
    coverages = [None] * logo_length

    for pos in range(logo_length):
        col_bases = _extract_column(seqs, pos)
        information_dict, coverage = _calc_frequencies(col_bases, base_vocabulary)

        coverages[pos] = coverage

        for base in base_vocabulary:
            matrix_columns[base][pos] = information_dict[base]
        # end for
    # end for

    final_matrix = pd.DataFrame(matrix_columns)
    final_matrix.index.rename('pos')

    return final_matrix, coverages
# end def create_matrix
