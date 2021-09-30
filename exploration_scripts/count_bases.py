#!/usr/bin/env python3
# -*- encoding: utf-8 -*-


from Bio import SeqIO

infpath = '/mnt/1.5_drive_0/16S_scrubbling/bacteria/gene_seqs/bacteria_pure_gene_seqs_annotated.fasta'
outfpath = '/mnt/1.5_drive_0/16S_scrubbling/bases_count.tsv'


with open(infpath, 'rt') as infile, open(outfpath, 'w') as outfile:

    outfile.write(f'seqID\ta\tt\tg\tc\tr\ty\tw\ts\tk\tm\th\tv\tb\td\tn\tlen\n')

    for i, record in enumerate(SeqIO.parse(infile, 'fasta')):

        seq = str(record.seq).upper()

        a = seq.count('A')
        t = seq.count('T')
        g = seq.count('G')
        c = seq.count('C')
        r = seq.count('R')
        y = seq.count('Y')
        w = seq.count('W')
        s = seq.count('S')
        k = seq.count('K')
        m = seq.count('M')
        h = seq.count('H')
        v = seq.count('V')
        b = seq.count('B')
        d = seq.count('D')
        n = seq.count('N')

        outfile.write(f'{record.id}\t{a}\t{t}\t{g}\t{c}\t{r}\t{y}\t{w}\t{s}\t{k}\t{m}\t{h}\t{v}\t{b}\t{d}\t{n}\t{len(seq)}\n')
    # end for
# end with

print('\nCompleted!')
print(outfpath)
