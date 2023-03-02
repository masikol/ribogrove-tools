
import sys
from io import StringIO
import subprocess as sp

import pandas as pd

from src.formatting import format_int_number


def _init_size_dict():

    ribogrove_size_dict = {
        'gene_num': {
            'Bacteria': None,
            'Archaea':  None,
            'Total':    None,
        },
        'uniq_gene_num': {
            'Bacteria': None,
            'Archaea':  None,
            'Total':    None,
        },
        'species_num': {
            'Bacteria': None,
            'Archaea':  None,
            'Total':    None,
        },
        'genome_num': {
            'Bacteria': None,
            'Archaea':  None,
            'Total':    None,
        },
        'cat1_genome_num': {
            'Bacteria': None,
            'Archaea':  None,
            'Total':    None,
        },
        'cat2_genome_num': {
            'Bacteria': None,
            'Archaea':  None,
            'Total':    None,
        },
        'cat3_genome_num': {
            'Bacteria': None,
            'Archaea':  None,
            'Total':    None,
        },
    }

    return ribogrove_size_dict
# end def


def count_dedup_seqs(fasta_fpath, seqkit_fpath, filter_str=None):

    if not filter_str is None:
        cmd = [
            seqkit_fpath, 'grep', '-nrp', f'";d__{filter_str};"', fasta_fpath,
            '|',
            seqkit_fpath, 'rmdup', '-s', 
            '|',
            seqkit_fpath, 'stats', '-T',
        ]
    else:
        cmd = [
            seqkit_fpath, 'rmdup', '-s', fasta_fpath,
            '|',
            seqkit_fpath, 'stats', '-T',
        ]
    # end if

    cmd_str = ' '.join(cmd)
    print(cmd_str)

    pipe = sp.Popen(cmd_str, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error while running command `{}`'.format(cmd_str))
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(1)
    # end if

    stdout_handle = StringIO(stdout_stderr[0].decode('utf-8'))
    tmp_stats_df = pd.read_csv(stdout_handle, sep='\t')

    print(tmp_stats_df)

    seq_count = int(tmp_stats_df.loc[0, 'num_seqs'])

    return seq_count
# end def


def make_ribogrove_size_dict(final_fasta_fpath, gene_stats_df, seqkit_fpath):
    ribogrove_size_dict = _init_size_dict()

    # Count Number of gene sequences
    ribogrove_size_dict['gene_num']['Bacteria'] = gene_stats_df[
        gene_stats_df['Domain'] == 'Bacteria'
    ]['seqID'].nunique()
    ribogrove_size_dict['gene_num']['Archaea'] = gene_stats_df[
        gene_stats_df['Domain'] == 'Archaea'
    ]['seqID'].nunique()
    ribogrove_size_dict['gene_num']['Total'] = gene_stats_df.shape[0]

    # Count Number of unique gene sequences
    ribogrove_size_dict['uniq_gene_num']['Bacteria'] = count_dedup_seqs(
        final_fasta_fpath,
        seqkit_fpath,
        'Bacteria'
    )
    ribogrove_size_dict['uniq_gene_num']['Archaea'] = count_dedup_seqs(
        final_fasta_fpath,
        seqkit_fpath,
        'Archaea'
    )
    ribogrove_size_dict['uniq_gene_num']['Total'] = count_dedup_seqs(
        final_fasta_fpath,
        seqkit_fpath
    )

    # Count Number of species
    ribogrove_size_dict['species_num']['Bacteria'] = gene_stats_df[
        gene_stats_df['Domain'] == 'Bacteria'
    ]['Species'].nunique()
    ribogrove_size_dict['species_num']['Archaea'] = gene_stats_df[
        gene_stats_df['Domain'] == 'Archaea'
    ]['Species'].nunique()
    ribogrove_size_dict['species_num']['Total'] = gene_stats_df['Species'].nunique()

    # Count Number of genomes
    ribogrove_size_dict['genome_num']['Bacteria'] = gene_stats_df[
        gene_stats_df['Domain'] == 'Bacteria'
    ]['asm_acc'].nunique()
    ribogrove_size_dict['genome_num']['Archaea'] = gene_stats_df[
        gene_stats_df['Domain'] == 'Archaea'
    ]['asm_acc'].nunique()
    ribogrove_size_dict['genome_num']['Total'] = gene_stats_df['asm_acc'].nunique()

    # Count Number of genomes of category 1,2,3
    for category in (1, 2, 3):
        ribogrove_size_dict[f'cat{category}_genome_num']['Bacteria'] = gene_stats_df[
            (gene_stats_df['Domain'] == 'Bacteria') & (gene_stats_df['category'] == category)
        ]['asm_acc'].nunique()
        ribogrove_size_dict[f'cat{category}_genome_num']['Archaea'] = gene_stats_df[
            (gene_stats_df['Domain'] == 'Archaea') & (gene_stats_df['category'] == category)
        ]['asm_acc'].nunique()
        ribogrove_size_dict[f'cat{category}_genome_num']['Total'] = gene_stats_df[
            gene_stats_df['category'] == category
        ]['asm_acc'].nunique()
    # end for

    for k, v in ribogrove_size_dict.items():
        print(f'{k}: {v}')
    # end for

    return ribogrove_size_dict
# end def


def format_size_dict(ribogrove_size_dict, thousand_separator, decimal_separator):

    fmt_ribogrove_size_dict = dict()

    for metric_name in ribogrove_size_dict.keys():
        fmt_ribogrove_size_dict[metric_name] = dict()
        for column_name, number in ribogrove_size_dict[metric_name].items():
            fmt_ribogrove_size_dict[metric_name][column_name] = format_int_number(
                number,
                thousand_separator
            )
        # end for
    # end for

    return fmt_ribogrove_size_dict
# end def
