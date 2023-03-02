#!/usr/bin/env python3

# Bacteria;Proteobacteria;Gammaproteobacteria;Pseudomonadales;Pseudomonadaceae;Pseudomonas;
# GCF_000287315.1:NC_018418.1:108826-110361:minus ;d__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Oceanospirillales;f__Halomonadaceae;g__Candidatus_Carsonella;s__ruddii; category:2

import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import re
import sys
import gzip
import argparse

from Bio import SeqIO


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-i',
    '--in-ribogrove-fasta',
    help='input RiboGrove-formatted fasta file (may be gzipped)',
    required=True
)

# Parameters

parser.add_argument(
    '-w',
    '--with-species',
    help='include species in the taxonomy',
    required=False,
    action='store_true'
)

# Output files

parser.add_argument(
    '-o',
    '--out-dada2-fasta',
    help='input dada2-train-set-formatted fasta.gz file',
    required=True
)

args = parser.parse_args()


# For convenience
infpath = os.path.abspath(args.in_ribogrove_fasta)
with_species = not args.with_species is None and args.with_species != False
outfpath = os.path.abspath(args.out_dada2_fasta)


DOMAIN_PATTERN  = re.compile(r';d__([A-Z][^;]+);')
PHYLUM_PATTERN  = re.compile(r';p__([A-Z][^;]+);')
CLASS_PATTERN   = re.compile(r';c__([A-Z][^;]+);')
ORDER_PATTERN   = re.compile(r';o__([A-Z][^;]+);')
FAMILY_PATTERN  = re.compile(r';f__([A-Z][^;]+);')
GENUS_PATTERN   = re.compile(r';g__([A-Z][^;]+);')
SPECIES_PATTERN = re.compile(r';s__([^;]+);')



def make_dada2_header_with_species(record):
    descr = record.description

    domain_name  = re.search(DOMAIN_PATTERN,  descr).group(1)
    phylum_name  = re.search(PHYLUM_PATTERN,  descr).group(1)
    class_name   = re.search(CLASS_PATTERN,   descr).group(1)
    order_name   = re.search(ORDER_PATTERN,   descr).group(1)
    family_name  = re.search(FAMILY_PATTERN,  descr).group(1)
    genus_name   = re.search(GENUS_PATTERN,   descr).group(1)
    species_name = re.search(SPECIES_PATTERN, descr).group(1)

    if phylum_name == 'NA':
        header = '{};'.format(domain_name)
    elif class_name == 'NA':
        header = '{};{};'.format(domain_name, phylum_name)
    elif order_name == 'NA':
        header = '{};{};{};'.format(domain_name, phylum_name, class_name)
    elif family_name == 'NA':
        header = '{};{};{};{};'.format(domain_name, phylum_name, class_name, order_name)
    elif genus_name == 'NA':
        header = '{};{};{};{};{};'.format(domain_name, phylum_name, class_name, order_name, family_name)
    elif species_name == 'NA':
        header = '{};{};{};{};{};{};'.format(domain_name, phylum_name, class_name, order_name, family_name, genus_name)
    else:
        header = '{};{};{};{};{};{};{};'.format(domain_name, phylum_name, class_name, order_name, family_name, genus_name, species_name)
    # end if

    header = header.replace('_', ' ')

    return header
# end def


def make_dada2_header_no_species(record):
    descr = record.description

    domain_name  = re.search(DOMAIN_PATTERN,  descr).group(1)
    phylum_name  = re.search(PHYLUM_PATTERN,  descr).group(1)
    class_name   = re.search(CLASS_PATTERN,   descr).group(1)
    order_name   = re.search(ORDER_PATTERN,   descr).group(1)
    family_name  = re.search(FAMILY_PATTERN,  descr).group(1)
    genus_name   = re.search(GENUS_PATTERN,   descr).group(1)

    if phylum_name == 'NA':
        header = '{};'.format(domain_name)
    elif class_name == 'NA':
        header = '{};{};'.format(domain_name, phylum_name)
    elif order_name == 'NA':
        header = '{};{};{};'.format(domain_name, phylum_name, class_name)
    elif family_name == 'NA':
        header = '{};{};{};{};'.format(domain_name, phylum_name, class_name, order_name)
    elif genus_name == 'NA':
        header = '{};{};{};{};{};'.format(domain_name, phylum_name, class_name, order_name, family_name)
    else:
        header = '{};{};{};{};{};{};'.format(domain_name, phylum_name, class_name, order_name, family_name, genus_name)
    # end if

    header = header.replace('_', ' ')

    return header
# end def


# Select formatting function
if with_species:
    make_dada2_header = make_dada2_header_with_species
    print('INFO: species names will be included')
else:
    make_dada2_header = make_dada2_header_no_species
    print('INFO: species names will not be included')
# end if


# == Proceed ==

if infpath.endswith('.gz'):
    open_input_fasta = gzip.open
else:
    open_input_fasta = open
# end if

print('Counting sequences...')
with open_input_fasta(infpath, 'rt') as infile:
    num_seqs = len(
        tuple(SeqIO.parse(infile, 'fasta'))
    )
# end with
print('{:,} sequences'.format(num_seqs))


print('Reformatting...')

step_frac = 0.05
step_num = round(num_seqs * step_frac)
next_update_num = step_num

sys.stdout.write('0/{} (0.0%)'.format(num_seqs))
sys.stdout.flush()


with open_input_fasta(infpath, 'rt') as infile, \
     gzip.open(outfpath, 'wt') as outfile:

    in_seq_records = SeqIO.parse(infile, 'fasta')

    for i, record in enumerate(in_seq_records):

        header = make_dada2_header(record)
        outfile.write(
            '>{}\n{}\n'.format(header, record.seq)
        )

        if i + 1 == next_update_num:
            percent = round(
                (i+1) / num_seqs * 100,
                2
            )
            sys.stdout.write(
                '\r{}/{} ({}%)'.format(i+1, num_seqs, percent)
            )
            next_update_num += step_num
            sys.stdout.flush()
        # end if
    # end for
# end with

print('\r{}/{} (100.0%)'.format(num_seqs, num_seqs))

print('\nOutput file: `{}`'.format(outfpath))
print('Completed!')
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')