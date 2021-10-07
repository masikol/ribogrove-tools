#!/usr/bin/env python3

import os
import sys
import math
import operator
from functools import reduce

import numpy as np
import pandas as pd
from Bio import SeqIO


def _check_lengths(seq_records) -> int:
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
    return tuple(
        map(
            lambda x: x[i],
            seqs
        )
    )
# end def _extract_column


def _calc_information(col_bases, base_vocabulary):

    aln_column = list(
        ''.join(col_bases).replace('-', '')
    )

    n_seqs = len(aln_column)

    max_information = math.log(len(base_vocabulary), 2)

    freq_dict = {base: aln_column.count(base) / n_seqs for base in base_vocabulary}

    # Calculate information
    # abs instead of minus in order not to allow "-0.0" values
    information = max_information - abs(
        reduce(
            operator.add,
            (
                freq * math.log(freq, 2) \
                    for base, freq in filter(
                        lambda x: x[0] in aln_column,
                        freq_dict.items()
                    )
            )
        )
    )

    information_dict = {base: 0.0 for base in base_vocabulary}

    # print(aln_column)
    # print(f'information = {information}')

    for base in information_dict.keys():
        base_freq = freq_dict[base]
        # print(f'base_freq({base}) = {base_freq}')
        if base_freq > 1e-9:
            # print(base_freq)
            information_dict[base] = information - ((-1) * base_freq * math.log(base_freq, 2))
        # end if
    # end for

    return information_dict
# end def _calc_information


def create_matrix(fasta_fpath):

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
    # matrix_columns['pos'] = np.array(range(logo_length))

    # print(matrix_columns)

    for pos in range(logo_length):
        col_bases = _extract_column(seqs, pos)
        information_dict = _calc_information(col_bases, base_vocabulary)

        for base in base_vocabulary:
            matrix_columns[base][pos] = information_dict[base]
        # end for

        # print(information_dict)
        # print('-----')
    # end for

    final_matrix = pd.DataFrame(matrix_columns)
    final_matrix.index.rename('pos')
    # print(final_matrix)

    # print('ok')

    return final_matrix
# end def create_matrix

