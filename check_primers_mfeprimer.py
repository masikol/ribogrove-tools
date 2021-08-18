#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

import os
import glob
import json

from Bio import SeqIO
import pandas as pd


primers = {
    '27F': 'AGAGTTTGATCMTGGCTCAG',
    '16S_BV2f': 'AGTGGCGGACGGGTGAGTAA',
    'CCF': 'CCAGACTCCTACGGGAGGCAGC',
    '16SAmpliconPCRForwardPrimer': 'CCTACGGGNGGCWGCAG',
    '336F-seq': 'ACGGYCCAGACTCCTACG',
    '338R': 'GCTGCCTCCCGTAGGAGT',
    '515F': 'GTGCCAGCMGCCGCGGTAA',
    '522R-seq': 'TATTACCGCGGCTGCTGGCAC',
    '16S_BV3r': 'CCGCGGCTGCTGGCAC',
    '519R': 'GWATTACCGCGGCKGCTG',
    '806R': 'GGACTACHVGGGTWTCTAAT',
    '784F': 'AGGATTAGATACCCTGGTA',
    '16SAmpliconPCRReversePrimer': 'GACTACHVGGGTATCTAATCC',
    '907R': 'CCGYCAATTCMTTTRAGTTT',
    '918F-seq': 'ACTCAAAKGAATTGACGGG',
    '926R-seq': 'CCGTCAATTCATTTGAGTTT',
    '968F': 'AACGCGAAGAACCTTAC',
    'CDR': 'CTTGTGCGGGCCCCCGTCAATTC',
    '1061R': 'CRRCACGAGCTGACGAC',
    '1100R': 'AGGGTTGCGCTCGTTG',
    '1114F': 'GCAACGAGCGCAACCC',
    '1193R': 'ACGTCATCCCCACCTTCC',
    '1237F': 'GGGCTACACACGYGCWAC',
    '1391R': 'GACGGGCGGTGTGTRCA',
    '1492Rl': 'GGTTACCTTGTTACGACTT',
}

primer_pairs = [
    # Full
    ['27F', '1492Rl'],
    # Two
    # V1-V2
    ['27F', '338R'],
    # V2-V3
    ['16S_BV2f', '16S_BV3r'],
    ['16S_BV2f', '519R'],
    ['16S_BV2f', '522R-seq'],
    # V3-V4
    ['16SAmpliconPCRForwardPrimer', '16SAmpliconPCRReversePrimer'],
    ['16SAmpliconPCRForwardPrimer', '806R'],
    ['CCF', '16SAmpliconPCRReversePrimer'],
    ['CCF', '806R'],
    ['336F-seq', '16SAmpliconPCRReversePrimer'],
    ['336F-seq', '806R'],
    # V4-V5
    ['515F', '907R'],
    ['515F', 'CDR'],
    ['515F', '926R-seq'],
    # V5-V6
    ['784F', '1061R'],
    ['784F', '1100R'],
    # V6-V7
    ['968F', '1193R'],
    ['918F-seq', '1193R'],
    # V7-V8
    ['1114F', '1391R'],
    # Three
    # V1-V3
    ['27F', '16S_BV3r'],
    ['27F', '519R'],
    ['27F', '522R-seq'],
    # V2-V4
    ['16S_BV2f', '16SAmpliconPCRReversePrimer'],
    ['16S_BV2f', '806R'],
    # V3-V5
    ['16SAmpliconPCRForwardPrimer', '907R'],
    ['16SAmpliconPCRForwardPrimer', 'CDR'],
    ['CCF', '907R'],
    ['CCF', 'CDR'],
    ['336F-seq', '907R'],
    ['336F-seq', 'CDR'],
    ['16SAmpliconPCRForwardPrimer', '926R-seq'],
    ['CCF', '926R-seq'],
    # V4-V6
    ['515F', '1061R'],
    ['515F', '1100R'],
    # V5-V7
    ['784F', '1193R'],
    # V6-V8
    ['968F', '1391R'],
    ['918F-seq', '1391R'],
    # V7-V9
    ['1114F', '1492Rl'],
]


mfeprimer = '/home/cager/Misc_soft/mfeprimer'
fasta_fpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/pure_genes_seqs_annotated.fasta'
categories_fpath = '/mnt/1.5_drive_0/16S_scrubbling/categories/bacteria_16S_genes_categories.tsv'

outdpath = '/mnt/1.5_drive_0/16S_scrubbling/check_primers/mfeprimer_results'

tmp_dir = '/mnt/1.5_drive_0/16S_scrubbling/check_primers/mfeprimer_results/tmp'
tmp_fasta = os.path.join(tmp_dir, 'tmpQ.fasta')
tmp_primers = os.path.join(tmp_dir, 'tmpP.fasta')
tmp_out_base = os.path.join(tmp_dir, 'tmpOUT')
tmp_out_json = f'{tmp_out_base}.json'
k_mer_size = 8

for d in (outdpath, tmp_dir):
    if not os.path.exists(d):
        os.makedirs(d)
    # end if
# end for

out_colnames = (
    'ass_id', 'seqID', 'product_size', 'ppc',
    'f_size', 'f_start', 'f_end', 'f_tm', 'f_dg', 'f_bind_len', 'f_pident', 'f_cover',
    'r_size', 'r_start', 'r_end', 'r_tm', 'r_dg', 'r_bind_len', 'r_pident', 'r_cover',
)


seqID_2_assID_dict = {
    row['seqID']: row['ass_id'] \
    for _, row in pd.read_csv(categories_fpath, sep='\t').iterrows()
}


for nameF, nameR in primer_pairs:
    outfpath = os.path.join(outdpath, f'{nameF}-{nameR}.tsv')
    with open(outfpath, 'w') as outfile:
        outfile.write('{}\n'.format('\t'.join(out_colnames)))
    # end with
# end for

num_seqs = len(tuple(SeqIO.parse(fasta_fpath, 'fasta')))
seq_records = tuple(SeqIO.parse(fasta_fpath, 'fasta'))

# accs = ['NZ_CP011342.2:c3404914-3403361',]


for i, seq_record in enumerate(seq_records):

    print(f'\rDoing {i+1}/{num_seqs}: {seq_record.id}', end=' ')

    ass_id = seqID_2_assID_dict[seq_record.id]

    for f in glob.iglob(f'{tmp_dir}/*'):
        os.unlink(f)
    # end for

    with open(tmp_fasta, 'wt') as tmp_fasta_file:
        tmp_fasta_file.write(f'>{seq_record.id}\n{seq_record.seq}\n')
    # end with

    index_cmd = f'{mfeprimer} index -i {tmp_fasta} -k {k_mer_size} -c 6'
    os.system(index_cmd)

    for nameF, nameR in primer_pairs:

        seqF = primers[nameF]
        seqR = primers[nameR]

        outfpath = os.path.join(outdpath, f'{nameF}-{nameR}.tsv')

        with open(tmp_primers, 'w') as tmp_file:
            tmp_file.write(f'>{nameF}\n{seqF}\n>{nameR}\n{seqR}\n')
        # end with

        cmd = f'{mfeprimer} spec --misEnd {k_mer_size} -k {k_mer_size} -c 6 -i {tmp_primers} -d {tmp_fasta} -j -o {tmp_out_base} 2> /dev/null'
        os.system(cmd)

        mfe_json = json.loads(open(tmp_out_json, 'rt').read())
        amp_list = mfe_json['AmpList']

        if not amp_list is None:
            # amp_count = len(amp_list)
            with open(outfpath, 'at') as outfile:

                for amp in amp_list:

                    product_size = amp['P']['Size']
                    ppc = amp['PPC']

                    f_size = amp['F']['Size']
                    f_start = amp['F']['Start']
                    f_end = amp['F']['End']
                    f_tm = amp['F']['Tm']
                    f_dg = amp['F']['Dg']
                    f_bind_len = len(amp['F']['Sseq'])
                    f_pident = amp['F']['Aseq'].count(':') / len(amp['F']['Aseq'])
                    f_cover = f_bind_len / f_size

                    r_size = amp['R']['Size']
                    r_start = amp['R']['Start']
                    r_end = amp['R']['End']
                    r_tm = amp['R']['Tm']
                    r_dg = amp['R']['Dg']
                    r_bind_len = len(amp['R']['Sseq'])
                    r_pident = amp['R']['Aseq'].count(':') / len(amp['R']['Aseq'])
                    r_cover = r_bind_len / r_size

                    out_tuple = map(
                        str,
                        (
                            ass_id, seq_record.id, product_size, ppc,
                            f_size, f_start, f_end, f_tm, f_dg, f_bind_len, f_pident, f_cover,
                            r_size, r_start, r_end, r_tm, r_dg, r_bind_len, r_pident, r_cover,
                        )
                    )

                    outfile.write('{}\n'.format('\t'.join(out_tuple)))
                # end for
            # end with
        # end if

        # os.unlink(tmp_primers)
        os.unlink(tmp_out_base)
        os.unlink(tmp_out_json)
    # end for
# end for

print('\nCompleted!')
print(outdpath)
