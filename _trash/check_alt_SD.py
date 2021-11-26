#!/usr/bin/env python3

import os
import sys
import gzip

import pandas as pd
from Bio import SeqIO
from Bio.SeqFeature import CompoundLocation, SeqFeature, FeatureLocation

# Bacteria
# ass_id_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/antiSD/alt_SD_check/no_antiSD_ass_IDs.txt'
# ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/bacteria_ass_refseq_accs_merged.tsv'
# tax_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/taxonomy/bacteria_per_genome_taxonomy.tsv'
# gbk_dpath = '/mnt/1.5_drive_0/16S_scrubbling/genomes-data/gbk'

# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/antiSD/alt_SD_check/upstream_seqs_CDSs.tsv'

# Archaea
# ass_id_fpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/antiSD/alt_SD_check/no_antiSD_ass_IDs.txt'
# ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/archaea_refseq_accs_merged.tsv'
# tax_fpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/taxonomy/archaea_per_genome_taxonomy.tsv'
# gbk_dpath = '/mnt/1.5_drive_0/16S_scrubbling/genomes-data/gbk'

# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/archaea/antiSD/alt_SD_check/upstream_seqs_CDSs.tsv'

# Test
# ass_id_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/antiSD/alt_SD_check/TEST_no_antiSD_ass_IDs.txt'
# ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/bacteria_ass_refseq_accs_merged.tsv'
# tax_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/taxonomy/bacteria_per_genome_taxonomy.tsv'
# gbk_dpath = '/mnt/1.5_drive_0/16S_scrubbling/genomes-data/gbk'

# outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/antiSD/alt_SD_check/TEST_upstream_seqs_CDSs.tsv'

# Test Bacillus subtilis (BACTEST)
ass_id_fpath = '/mnt/1.5_drive_0/16S_scrubbling/antiSD/bacteria/alt_SD_check/BACTEST_no_antiSD_ass_IDs.txt'
ass_acc_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/bacteria_ass_refseq_accs_merged.tsv'
tax_fpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/taxonomy/bacteria_per_genome_taxonomy.tsv'
gbk_dpath = '/mnt/1.5_drive_0/16S_scrubbling/genomes-data/gbk'

outfpath = '/mnt/1.5_drive_0/16S_scrubbling/antiSD/bacteria/alt_SD_check/BACTEST_upstream_seqs_CDSs.tsv'


max_spacer_size = 15
sd_len = 6
upstream_seq_size = max_spacer_size + sd_len

with open(ass_id_fpath, 'rt') as ass_id_file:
    all_ass_ids = set(
        map(
            lambda x: int(x.strip()),
            ass_id_file.readlines()
        )
    )
# end with

ass_acc_df = pd.read_csv(ass_acc_fpath, sep='\t').query('ass_id in @all_ass_ids')
# print(ass_acc_df.shape)
# print(ass_acc_df.head())

tax_df = pd.read_csv(tax_fpath, sep='\t').query('ass_id in @all_ass_ids')
# print(tax_df.shape)
# print(tax_df.head())

no_antiSD_df = ass_acc_df.merge(
    tax_df.drop(['accs', 'taxID', 'tax_name'], axis=1),
    on='ass_id',
    how='left'
)


# print(no_antiSD_df.shape)
# print(no_antiSD_df.head())

ass_ids = set(no_antiSD_df['ass_id'])
# ass_ids = {
#     1020931,
#     500738,
#     6720021,
#     6720041,
#     31868,
#     1865971,
#     421728,
#     1023351,
#     1166741,
#     8103341,
#     10079181,
# }

with open(outfpath, 'wt') as outfile:

    outfile.write('ass_id\tacc\tlocus_tag\tstart\tend\tupstream_seq\tjoint_loc\n')

    for i, ass_id in enumerate(ass_ids):

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        curr_no_antiSD_df = no_antiSD_df[no_antiSD_df['ass_id'] == ass_id]
        accs = tuple(curr_no_antiSD_df['acc'])

        # print(accs)

        for acc in accs:

            gbk_fpath = os.path.join(gbk_dpath, f'{acc}.gbk.gz')

            with gzip.open(gbk_fpath, 'rt') as gbk_file:
                seq_record = tuple(SeqIO.parse(gbk_file, 'gb'))[0]
                # print(seq_record.id, seq_record.seq[:10])
            # end with

            for f in seq_record.features:
                if f.type == 'CDS':

                    # print(f'f.strand = {f.strand}')
                    # print(type(f.location))

                    # If location if "join"ed
                    joint_loc = 0
                    if isinstance(f.location, CompoundLocation):
                        joint_loc = 1
                        # print()
                        # print(f.location)
                        # print(f.location.parts)
                        if f.strand == 1:
                            location = f.location.parts[0]
                        else:
                            location = f.location.parts[1]
                        # end if
                    else:
                        location = f.location
                    # end if

                    # print(location)
                    replicon_len = len(seq_record.seq)

                    cross_start_feature = False

                    if f.strand == 1:
                        start = location.start
                        upstream_start = start - upstream_seq_size
                        upstream_end = start
                        if upstream_start >= replicon_len:
                            upstream_start -= replicon_len
                            cross_start_feature = True
                        # end if
                        if upstream_end >= replicon_len:
                            upstream_end -= replicon_len
                            cross_start_feature = True
                        # end if
                        # print(f'start = {start}')
                        if upstream_start >= 0:
                        # if upstream_start <= upstream_end:
                            # upstream_seq = seq_record.seq[upstream_start : upstream_end]
                            if upstream_start > upstream_end:
                                upstream_seq = SeqFeature(FeatureLocation(0, upstream_end)).extract(seq_record) \
                                             + SeqFeature(FeatureLocation(upstream_start, replicon_len)).extract(seq_record)
                            else:
                                upstream_seq = SeqFeature(FeatureLocation(upstream_start, upstream_end)).extract(seq_record)
                            # end if
                        else:
                            cross_start_feature = True
                            upstream_start = replicon_len - abs(upstream_start)
                            upstream_seq = SeqFeature(FeatureLocation(upstream_start, replicon_len)).extract(seq_record) \
                                         + SeqFeature(FeatureLocation(0, upstream_end)).extract(seq_record)
                        # end if
                    else:
                        end = location.end
                        upstream_start = end
                        upstream_end = end + upstream_seq_size
                        if upstream_start >= replicon_len:
                            upstream_start -= replicon_len
                            cross_start_feature = True
                        # end if
                        if upstream_end >= replicon_len:
                            cross_start_feature = True
                            upstream_end -= replicon_len
                        # end if
                        # print(f'end = {end}')
                        if upstream_start >= 0:
                        # if upstream_start <= upstream_end:
                            # upstream_seq = seq_record.seq[upstream_start : upstream_end]
                            if upstream_start > upstream_end:
                                upstream_seq = SeqFeature(FeatureLocation(0, upstream_end)).extract(seq_record).reverse_complement() \
                                             + SeqFeature(FeatureLocation(upstream_start, replicon_len)).extract(seq_record).reverse_complement()
                            else:
                                upstream_seq = SeqFeature(FeatureLocation(upstream_start, upstream_end)).extract(seq_record).reverse_complement()
                            # end if
                        else:
                            cross_start_feature = True
                            upstream_start = replicon_len - abs(upstream_start)
                            upstream_seq = SeqFeature(FeatureLocation(upstream_start, replicon_len)).extract(seq_record).reverse_complement() \
                                         + SeqFeature(FeatureLocation(0, upstream_end)).extract(seq_record).reverse_complement()
                        # end if
                        # upstream_seq = seq_record.seq[upstream_start : upstream_end].reverse_complement()
                        # upstream_seq = SeqFeature(FeatureLocation(upstream_start, upstream_end)).extract(seq_record).reverse_complement()
                    # end if

                    if not (cross_start_feature and seq_record.annotations['topology'].lower() != 'circular'):
                        locus_tag = f.qualifiers['locus_tag'][0]
                        upstream_seq = str(upstream_seq.seq)

                        out_start = str(upstream_start+1).replace('>', '').replace('<', '')
                        out_end = str(upstream_end).replace('>', '').replace('<', '')
                        if out_end == '0':
                            out_end = str(len(seq_record.seq))
                        # end if

                        outfile.write(f'{ass_id}\t{acc}\t{locus_tag}\t{out_start}\t{out_end}\t{upstream_seq}\t{joint_loc}\n')
                    # end if

                    # print(f)
                    # print(upstream_seq)
                    # input()
                # end if
            # end for
        # end for

    # end for
# ednd with

print('\nCompleted!')
print(outfpath)
