import sys
from io import StringIO
import subprocess as sp

import polars as pl

from src.formatting import format_int_number
from src.util import is_validlike_species_name


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


def _count_unique(gene_stats_df, column):
    return gene_stats_df.select(
        pl.col(column).drop_nulls().n_unique()
    ).item()
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
    tmp_stats_df = pl.read_csv(stdout_handle, separator='\t')

    print(tmp_stats_df)

    seq_count = int(tmp_stats_df.item(0, 'num_seqs'))

    return seq_count
# end def


def make_ribogrove_size_dict(final_fasta_fpath, gene_stats_df, seqkit_fpath):
    ribogrove_size_dict = _init_size_dict()

    bacteria_df = gene_stats_df.filter(pl.col('Domain') == 'Bacteria')
    archaea_df = gene_stats_df.filter(pl.col('Domain') == 'Archaea')

    # Count Number of gene sequences
    ribogrove_size_dict['gene_num']['Bacteria'] = _count_unique(
        bacteria_df, 'seqID'
    )
    ribogrove_size_dict['gene_num']['Archaea'] = _count_unique(
        archaea_df, 'seqID'
    )
    ribogrove_size_dict['gene_num']['Total'] = gene_stats_df.height

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
    total_species_count = 0
    for domain, domain_df in (('Bacteria', bacteria_df), ('Archaea', archaea_df),):
        species_set = frozenset(domain_df['Species'].to_list())
        ribogrove_size_dict['species_num'][domain] = len(tuple(filter(
            is_validlike_species_name,
            species_set
        )))
        total_species_count += ribogrove_size_dict['species_num'][domain]
    # end for
    ribogrove_size_dict['species_num']['Total'] = total_species_count

    # Count Number of genomes
    ribogrove_size_dict['genome_num']['Bacteria'] = _count_unique(
        bacteria_df, 'asm_acc'
    )
    ribogrove_size_dict['genome_num']['Archaea'] = _count_unique(
        archaea_df, 'asm_acc'
    )
    ribogrove_size_dict['genome_num']['Total'] = _count_unique(
        gene_stats_df, 'asm_acc'
    )

    # Count Number of genomes of category 1,2,3
    for category in (1, 2, 3):
        category_df = gene_stats_df.filter(pl.col('category') == category)
        ribogrove_size_dict[f'cat{category}_genome_num']['Bacteria'] = _count_unique(
            category_df.filter(pl.col('Domain') == 'Bacteria'), 'asm_acc'
        )
        ribogrove_size_dict[f'cat{category}_genome_num']['Archaea'] = _count_unique(
            category_df.filter(pl.col('Domain') == 'Archaea'), 'asm_acc'
        )
        ribogrove_size_dict[f'cat{category}_genome_num']['Total'] = _count_unique(
            category_df, 'asm_acc'
        )
    # end for

    for k, v in ribogrove_size_dict.items():
        print(f'{k}: {v}')
    # end for

    return ribogrove_size_dict
# end def


def format_size_dict(ribogrove_size_dict, thousand_separator, decimal_separator):

    fmt_ribogrove_size_dict = {}

    for metric_name in ribogrove_size_dict.keys():
        fmt_ribogrove_size_dict[metric_name] = {}
        for column_name, number in ribogrove_size_dict[metric_name].items():
            fmt_ribogrove_size_dict[metric_name][column_name] = format_int_number(
                number,
                thousand_separator
            )
        # end for
    # end for

    return fmt_ribogrove_size_dict
# end def
