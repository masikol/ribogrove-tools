#!/usr/bin/env python3
# -*- encoding: utf-8 -*-


from Bio import SeqIO

infpath = '/mnt/1.5_drive_0/16S_scrubbling/gene_seqs/pure_genes_seqs_annotated.fasta'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bases_count.tsv'


with open(infpath, 'rt') as infile, open(outfpath, 'w') as outfile:

    outfile.write(f'seqID\ta\tt\tg\tc\tlen\n')

    for i, record in enumerate(SeqIO.parse(infile, 'fasta')):

        seq = str(record.seq).upper()

        a = seq.count('A')
        t = seq.count('T')
        g = seq.count('G')
        c = seq.count('C')

        outfile.write(f'{record.id}\t{a}\t{t}\t{g}\t{c}\t{len(seq)}\n')
    # end for
# end with

print('\nCompleted!')
print(outfpath)
