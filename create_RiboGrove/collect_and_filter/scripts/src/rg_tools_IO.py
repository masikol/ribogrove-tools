
import sys
import gzip

import pandas as pd
import polars as pl
from Bio import SeqIO


def read_ass_sum_file(infpath, raw_summary=False):
    try:
        asm_sum_df = read_asm_sum_file_polars(infpath, raw_summary=raw_summary)
    except pl.exceptions.ComputeError as err:
        print('polars.read_csv error:', file=sys.stderr)
        print(err, file=sys.stderr)
        print('Trying to read using pandas...', file=sys.stderr)
        asm_sum_df = pl.from_pandas(
            read_asm_sum_file_pandas(infpath, raw_summary=raw_summary)
        )
        print('Success!', file=sys.stderr)
    # end try

    return asm_sum_df
# end def


def read_asm_sum_file_polars(infpath, raw_summary=False):

    columns = [
        '#assembly_accession',
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
        'asm_submitter',
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
        '#assembly_accession': pl.String,
        'bioproject': pl.String,
        'biosample': pl.String,
        'wgs_master': pl.String,
        'refseq_category': pl.String,
        'taxid': pl.UInt32,
        'species_taxid': pl.UInt32,
        'organism_name': pl.String,
        'infraspecific_name': pl.String,
        'isolate': pl.String,
        'version_status': pl.String,
        'assembly_level': pl.String,
        'release_type': pl.String,
        'genome_rep': pl.String,
        'seq_rel_date': pl.String,
        'asm_name': pl.String,
        'asm_submitter': pl.String,
        'gbrs_paired_asm': pl.String,
        'paired_asm_comp': pl.String,
        'ftp_path': pl.String,
        'excluded_from_refseq': pl.String,
        'relation_to_type_material': pl.String,
        'asm_not_live_date': pl.String,
        # 'assembly_type': pl.String,
        # 'group': pl.String,
        # 'genome_size': pl.String,
        # 'genome_size_ungapped': pl.String,
        # 'gc_percent': pl.String,
        # 'replicon_count': pl.String,
        # 'scaffold_count': pl.String,
        # 'contig_count': pl.String,
        # 'annotation_provider': pl.String,
        # 'annotation_name': pl.String,
        # 'annotation_date': pl.String,
        # 'total_gene_count': pl.String,
        # 'protein_coding_gene_count': pl.String,
        # 'non_coding_gene_count': pl.String,
        # 'pubmed_id': pl.String,
    }

    if not raw_summary:
        columns[0] = 'asm_acc'
        del colDTypes['#assembly_accession']
        colDTypes['asm_acc'] =  pl.String
    # end if

    asm_sum_df = pl.read_csv(
        infpath,
        separator='\t',
        comment_prefix='##',
        n_threads=1,
        has_header=True,
        columns=columns,
        schema_overrides=colDTypes,
        null_values=['NA', 'na', '']
    )

    if raw_summary:
        asm_sum_df = asm_sum_df.rename({
            '#assembly_accession': 'asm_acc',
        })
    # end if

    return asm_sum_df
# end def


def read_asm_sum_file_pandas(infpath, raw_summary=False):

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
        'asm_submitter',
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
        'asm_submitter': str,
        'gbrs_paired_asm': str,
        'paired_asm_comp': str,
        'ftp_path': str,
        'excluded_from_refseq': str,
        'relation_to_type_material': str,
        'asm_not_live_date': str,
        # 'assembly_type': pl.String,
        # 'group': pl.String,
        # 'genome_size': pl.String,
        # 'genome_size_ungapped': pl.String,
        # 'gc_percent': pl.String,
        # 'replicon_count': pl.String,
        # 'scaffold_count': pl.String,
        # 'contig_count': pl.String,
        # 'annotation_provider': pl.String,
        # 'annotation_name': pl.String,
        # 'annotation_date': pl.String,
        # 'total_gene_count': pl.String,
        # 'protein_coding_gene_count': pl.String,
        # 'non_coding_gene_count': pl.String,
        # 'pubmed_id': pl.String,
    }

    # If the summary file is raw (freshly downloaded),
    #   we need to skip the very first line of it.
    skiprows = [0] if raw_summary else 0

    with gzip.open(infpath, 'rt') as infile:
        asm_sum_df = pd.read_csv(
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
    return asm_sum_df
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
