
import sys
import gzip

import pandas as pd
from Bio import SeqIO


def read_ass_sum_file(infpath, raw_summary=False):

    colNames = [
        'asm_acc',
        'bioproject',
        'biosample',
        'wgs_master',
        'refseq_category',
        'taxid',
        'species_taxid',
        'organism_name',
        'infraspecific_name',
        'isolate',
        'version_status',
        'assembly_level',
        'release_type',
        'genome_rep',
        'seq_rel_date',
        'asm_name',
        'submitter',
        'gbrs_paired_asm',
        'paired_asm_comp',
        'ftp_path',
        'excluded_from_refseq',
        'relation_to_type_material',
        'asm_not_live_date',
        # 'assembly_type',
        # 'group',
        # 'genome_size',
        # 'genome_size_ungapped',
        # 'gc_percent',
        # 'replicon_count',
        # 'scaffold_count',
        # 'contig_count',
        # 'annotation_provider',
        # 'annotation_name',
        # 'annotation_date',
        # 'total_gene_count',
        # 'protein_coding_gene_count',
        # 'non_coding_gene_count',
        # 'pubmed_id',
    ]

    colDTypes = {
        'asm_acc': str,
        'bioproject': str,
        'biosample': str,
        'wgs_master': str,
        'refseq_category': str,
        'taxid': pd.UInt32Dtype(),
        'species_taxid': pd.UInt32Dtype(),
        'organism_name': str,
        'infraspecific_name': str,
        'isolate': str,
        'version_status': str,
        'assembly_level': str,
        'release_type': str,
        'genome_rep': str,
        'seq_rel_date': str,
        'asm_name': str,
        'submitter': str,
        'gbrs_paired_asm': str,
        'paired_asm_comp': str,
        'ftp_path': str,
        'excluded_from_refseq': str,
        'relation_to_type_material': str,
        'asm_not_live_date': str,
        # 'assembly_type': str,
        # 'group': str,
        # 'genome_size': str,
        # 'genome_size_ungapped': str,
        # 'gc_percent': str,
        # 'replicon_count': str,
        # 'scaffold_count': str,
        # 'contig_count': str,
        # 'annotation_provider': str,
        # 'annotation_name': str,
        # 'annotation_date': str,
        # 'total_gene_count': str,
        # 'protein_coding_gene_count': str,
        # 'non_coding_gene_count': str,
        # 'pubmed_id': str,
    }

    # If the summary file is raw (freshly downloaded),
    #   we need to skip the very first line of it.
    skiprows = [0] if raw_summary else 0

    with gzip.open(infpath, 'rt') as infile:
        ass_sum_df = pd.read_csv(
            infile,
            sep='\t',
            skiprows=skiprows,
            # engine='c',
            header=0,
            names=colNames,
            usecols=colNames,
            dtype=colDTypes,
            na_values=['NA', 'na', '']
        )
    # end with
    return ass_sum_df
# end def


def read_and_filter_fasta(in_fasta_fpath,
                          filter_fpaths=[],
                          blacklist=set(),
                          whitelist=set()):
    # This read_and_filter_fasta function reads a fasta file,
    #   removes sequences having ids listed in one of the filter files,
    #   removes sequences listed in the blacklist
    #   but keeps sequences listed in the whitelist.

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
