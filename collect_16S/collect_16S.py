#!/usr/bin/env python3

import os
import gzip
import statistics as sts

import pandas as pd
from Bio import SeqIO

# in_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/test_bacteria_ass_refseq_accs_merged.tsv'
in_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria_ass_refseq_accs_merged.tsv'
gbk_dpath = '/mnt/1.5_drive_0/preprocess-dev/own_db/bacteria/pileup/genomes-dwnld/genomes-data/gbk'

fasta_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected.fasta.gz'
outstats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/all_collected_collect_16S_stats.tsv'
# fasta_outfpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/test_all_collected.fasta.gz'
# outstats_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/test_all_collected_collect_16S_stats.tsv'

cmsearch = '/home/cager/Misc_soft/infernal/infernal-1.1.1/bin/cmsearch'
rfam_12_0_fpath = '/mnt/1.5_drive_0/16S_scrubbling/rfam/RF00177.12.0.cm'
tblout_header = 'target_name\taccession\tquery_name\taccession\tmdl\tmdl_from\tmdl_to\tseq_from\tseq_to\tstrand\ttrunc\tpass\tgc\tbias\tscore\tEvalue\tinc\tdescription_of_target'

# possible_seqstart_trunc_note = '16S ribosomal RNA rRNA prediction is too short'
# cmsearch_note = 'Derived by automated computational analysis using gene prediction method: cmsearch.'

# note:
# 16S ribosomal RNA rRNA prediction is too short
# note:
# possible 16S ribosomal RNA but does not have goodblast hits on one or both of the ends

stats_header = [
    'ass_id', 'refseq_id', 'acc', 'title',
    'seq_start_truncation', 'improper_16S_annotation', 'topology',
    'num_genes', 'min_len', 'max_len', 'mean_len', 'median_len',
]


ssu_product_names = {
    '16S RIBOSOMAL RNA',
    'SMALL SUBUNIT RIBOSOMAL RNA',
    'SMALL RIBOSOMAL RNA',
    'S-RRNA',
    'RIBOSOMAL RNA-16S'
}

trunc_ssu_notes = {
    '16S RIBOSOMAL RNA RRNA PREDICTION IS TOO SHORT',
    'POSSIBLE 16S RIBOSOMAL RNA',
}

def is_ssu(feature):
    qualifiers = feature.qualifiers

    norm_ssu = 'product' in qualifiers.keys() \
               and qualifiers['product'][0].upper() in ssu_product_names

    maybe_trunc_ssu = False
    if 'note' in qualifiers.keys():
        concat_notes = ''.join(qualifiers['note']).upper()
        maybe_trunc_ssu = any(map(lambda note: note in concat_notes, trunc_ssu_notes))
    # end if

    return norm_ssu or maybe_trunc_ssu
# end def is_ssu

def filter_ssu_genes(features):
    return tuple(
        filter(
            is_ssu,
            features
        )
    )
# end def filter_ssu_genes


def is_annotated_with_pgap(gbrecord):
    try:
        struct_comment = gbrecord.annotations['structured_comment']
    except KeyError as err:
        return False
    # end try

    if 'Genome-Annotation-Data' in struct_comment.keys():
        assembly_key = 'Genome-Annotation-Data'
    else:
        return False
    # end if

    assembly_data = struct_comment[assembly_key]

    if 'Annotation Pipeline' in assembly_data.keys():
        annot_pipe_key = 'Annotation Pipeline'
    else:
        return False
    # end if

    annot_pipe = assembly_data[annot_pipe_key]

    if 'NCBI PROKARYOTIC GENOME ANNOTATION PIPELINE' in annot_pipe.upper():
        return True
    else:
        return False
    # end if
# end def


# def need_to_run_cmsearch(features):

#     for f in features:
#         if not 'note' in f.qualifiers:
#             return True
#         elif not cmsearch_note in f.qualifiers['note']:
#             return True
#         # end if
#     # end for

#     return False
# # end def get_notes


# def remove_internal_truncated_genes(features, gbrecord):
#     seq_end_coord = len(gbrecord.seq)

#     def is_internal_truncated(feature):
#         truncated = 'note' in feature.qualifiers.keys() \
#                     and possible_seqstart_trunc_note in feature.qualifiers['note']
#         if not truncated:
#             return True
#         # end if
#         internal = feature.location.start != 0 and feature.location.end != seq_end_coord
#         return not (truncated and internal)
#     # end def is_internal_truncated

#     return list(
#         filter(
#             is_internal_truncated,
#             features
#         )
#     )
# # end def remove_internal_truncated_genes


def seq_start_may_truncate_ssu(gbrecord, features):
    for f in features:
        if f.location.start == 0 or f.location.end == len(gbrecord.seq):
            return True
        # end if
    # end for
    return False
# end def seq_start_may_truncate_ssu


def extract_gene_as_is(feature, gbrecord):

    seq_start = feature.location.start + 1
    seq_end = feature.location.end
    seq_strand = feature.location.strand

    seq = gbrecord.seq[seq_start-1 : seq_end]
    if seq_strand == -1:
        seq = seq.reverse_complement()
    # end if

    strand_str = 'plus' if seq_strand == 1 else 'minus'

    header = f'{gbrecord.id}:{seq_start}-{seq_end}_{strand_str} {gbrecord.description}'

    return header, seq
# end def extract_gene_as_is


def run_cmsearch(fasta_fpath):
    tblout_fpath = 'tmpXXX_tblout.tsv'
    out_fpath = 'tmpXXX_cmsearch_out.txt'
    # ~/Misc_soft/infernal/infernal-1.1.1/bin/cmsearch --noali -o cmsearch_output.txt \
    #   --tblout cmsearch_tblout.tsv --cpu 4 SSU_rRNA_bacteria_12.0.cm ../genomes-data/fasta/NZ_CP060094.1.fasta.gz
    cmd = f'{cmsearch} --noali -o {out_fpath} --tblout {tblout_fpath} --cpu 6 {rfam_12_0_fpath} {fasta_fpath}'
    # cmd = f'{cmsearch} -o {out_fpath} --tblout {tblout_fpath} --cpu 6 {rfam_12_0_fpath} {fasta_fpath}'

    exit_code = os.system(cmd)
    if exit_code != 0:
        raise OSError(f'exit code {exit_code}')
    # end if
    return tblout_fpath
# end def run_cmsearch


def reformat_tblout(tblout_fpath):

    with open(tblout_fpath, 'rt') as tblout_file:
        lines = list(
            map(
                str.strip,
                filter(
                    lambda x: x[0] != '#',
                    tblout_file.readlines()
                )
            )
        )
    # end with

    for i in range(len(lines)):
        for space_num in range(20, 1, -1):
            lines[i] = lines[i].replace(' '*space_num, ' ')
        # end for
    # end for

    for i in range(len(lines)):
        lines[i] = lines[i].replace(' ', '\t')
    # end for

    with open(tblout_fpath, 'wt') as tblout_file:
        tblout_file.write(f'{tblout_header}\n')
        tblout_file.write('\n'.join(lines) + '\n')
    # end with
# end def reformat_tblout


def extract_reannotated_genes(gbrecord, topology):

    original_len = len(gbrecord.seq)

    # Make circular sequence a bit longer in order to annotate properly
    #   genes truncated by sequecne start
    if topology == 'circular':
        # print(f'Len before = {original_len}')
        gbrecord.seq = gbrecord.seq + gbrecord.seq[:2000]
        # print(f'Len after = {len(gbrecord.seq)}')
    # end if

    tmp_fasta = 'tmpXXX.fasta'
    with open(tmp_fasta, 'w') as tmpf:
        tmpf.write(f'>{gbrecord.id}\n{str(gbrecord.seq)}\n')
    # end with
    tblout_fpath = run_cmsearch(tmp_fasta)
    reformat_tblout(tblout_fpath)

    tblout_df = pd.read_csv(tblout_fpath, sep='\t')
    # tblout_df = tblout_df[tblout_df['trunc'] == 'no']

    genes = list()

    for i, row in tblout_df.iterrows():
        seq_start = row['seq_from']
        seq_end = row['seq_to']
        seq_strand = row['strand']

        if seq_strand == '+':
            # seq_start -= 1
            seq = gbrecord.seq[seq_start : seq_end]
        else:
            seq_start, seq_end = seq_end, seq_start
            # seq_start -= 1
            seq = gbrecord.seq[seq_start : seq_end].reverse_complement()
        # end if

        strand_str = 'plus' if seq_strand == '+' else 'minus'

        # Amend coordinated of genes truncated by sequence start
        #   in order to replace coordinates greater than `original_len`
        seq_start_for_header = seq_start
        if seq_start_for_header > original_len:
            seq_start_for_header -= original_len
        # end if
        seq_end_for_header = seq_end
        if seq_end_for_header > original_len:
            seq_end_for_header -= original_len
        # end if

        seq_header = '{}:{}-{}_{} {}'\
            .format(
                gbrecord.id,
                seq_start_for_header,
                seq_end_for_header,
                strand_str,
                gbrecord.description
            )

        genes.append(
            [
                seq_header,
                seq
            ]
        )

    # end for

    # for g in genes:
    #     print(f'>{g[0]}')
    #     print(f'LEN = {len(g[1])}')

    return genes
# end def extract_reannotated_genes


def get_seq_lengths(extracted_genes):
    lengths = [len(gene[1]) for gene in extracted_genes]
    return lengths
# end def get_seq_lengths

def calc_gene_stats(extracted_genes):
    lengths = get_seq_lengths(extracted_genes)

    if len(lengths) == 0:
        return 0, 'NA', 'NA', 'NA', 'NA'
    # end if

    min_len = min(lengths)
    max_len = max(lengths)
    mean_len = float(sts.mean(lengths))
    median_len = float(sts.median(lengths))

    return len(lengths), min_len, max_len, mean_len, median_len
# end def calc_gene_stats



acc_df = pd.read_csv(
    in_acc_fpath,
    sep='\t'
)

n_accs = acc_df.shape[0]

with gzip.open(fasta_outfpath, 'wt') as fasta_outfile, open(outstats_fpath, 'wt') as stats_outfile:

    stats_outfile.write('{}\n'.format('\t'.join(stats_header)))

    for i, row in acc_df.iterrows():
        ass_id = row['ass_id']
        refseq_id = row['refseq_id']
        acc = row['acc']
        title = row['title']
        print(f'\rDoing {i+1}/{n_accs}: {acc}', end=' '*10)

        gbk_fpath = os.path.join(
            gbk_dpath,
            f'{acc}.gbk.gz'
        )

        with gzip.open(gbk_fpath, 'rt') as gbfile:
            gbrecord = tuple(SeqIO.parse(gbfile, 'gb'))[0]
        # end with

        # print()

        ssu_features = filter_ssu_genes(gbrecord.features)
        # ssu_features = remove_internal_truncated_genes(ssu_features, gbrecord)

        # for f in ssu_features:
        #     print(f.location.start, f.location.end)
        # # end for

        seq_start_truncation = seq_start_may_truncate_ssu(gbrecord, ssu_features)
        # print(f'seq_start_truncation = {seq_start_truncation}')

        improper_16S_annotation = not is_annotated_with_pgap(gbrecord)
        # improper_16S_annotation = need_to_run_cmsearch(ssu_features)
        # print(f'improper_16S_annotation = {improper_16S_annotation}')

        topology = None
        if 'topology' in gbrecord.annotations.keys():
            topology = gbrecord.annotations['topology']
        # end if

        # print(f'topology: {topology}')

        if improper_16S_annotation or (seq_start_truncation and topology == 'circular'):
            extracted_genes = extract_reannotated_genes(gbrecord, topology)
            for header, seq in extracted_genes:
                fasta_outfile.write(f'>{header}\n{seq}\n')
            # end for
        # if not improper_16S_annotation and not seq_start_truncation:
        else:
            extracted_genes = [extract_gene_as_is(f, gbrecord) for f in ssu_features]
            for header, seq in extracted_genes:
                fasta_outfile.write(f'>{header}\n{seq}\n')
            # end for
        # end if

        stats_outfile.write(f'{ass_id}\t{refseq_id}\t{acc}\t{title}\t')
        stats_outfile.write(f'{1 if seq_start_truncation else 0}\t{1 if improper_16S_annotation else 0}\t{topology}\t')

        num_genes, min_len, max_len, mean_len, median_len = calc_gene_stats(extracted_genes)
        stats_outfile.write(f'{num_genes}\t{min_len}\t{max_len}\t{mean_len}\t{median_len}\n')

        # print('-' * 20)
        # input()

        # fasta_outfile.write(f'{acc}\t{seqtech}\n')
    # end for
# end with

print('\nRunning seqkit rmdup...')
tmpfasta = os.path.join(
    os.path.dirname(fasta_outfpath),
    'tmpBILLY.fasta'
)
os.system(f'zcat {fasta_outfpath} | seqkit rmdup -n | gzip > {tmpfasta}')
os.system(f'zcat {tmpfasta} | seqkit seq -u | gzip > {fasta_outfpath}')
print('Done\n')

print('\nCompleted!')
print(fasta_outfpath)
print(outstats_fpath)
