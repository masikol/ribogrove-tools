#!/usr/bin/env python3

# The script extracts sequences of 16S genes from the downloaded genomes in GenBank format.
# If the GenBank file of a genomic sequence does not include the statement that it
#   was annotated with the PGAP (https://www.ncbi.nlm.nih.gov/genome/annotation_prok/),
#   the script will re-annotate this sequence with cmsearch (http://eddylab.org/infernal/) 1.1.1
#   and Rfam (https://rfam.xfam.org/) 12.0, to be fully consistent with the PGAP.
# If the original 16S rRNA gene annotation of a circular sequence includes the first or
#   the last base of the sequence, the script appends 2,000 5′‑terminal bases to 3′-end
#   of this sequence prior to re-annotation, in order not to miss genes interrupted by sequence start.

## Command line arguments
### Input files:
# 1. `-i / --assm-acc-file` -- a TSV file of 4 columns:(`ass_id`, `gi_number`, `acc`, `title`).
#   This file is the output of the script `merge_assID2acc_and_remove_WGS.py`. Mandatory.
# 2. `-g/--gbk-dir` -- the directory where the downloaded `.gbk.gz` files are located
#   (see script `download_genomes.py`). Mandatory.

### Output files:
# 1. `-o / --out-fasta` -- a fasta file containing sequences of extracted genes. Mandatory.
# 2. `-s / --out-stats` -- a file with per-replicon statistics of extracted 16S genes:
#   how many genes, minimum/maximum length etc. Mandatory.

### Dependencies:
# 1. `--cmsearch` -- a `cmsearch` program executable from Infernal (http://eddylab.org/infernal/).
#   Please, use Infernal 1.1.1 with this script, to be fully consistent with the PGAP. Mandatory.
# 2. `-r / --rfam-family-cm` -- an (uncompressed) `.cm` file containing a Rfam's (version 12.0)
#    covariance models: ftp://ftp.ebi.ac.uk/pub/databases/Rfam/12.0/Rfam.cm.gz. Mandatory.
# 3. `--seqkit` -- a `seqkit` executable: github.com/shenwei356/seqkit.
#   Mandatory.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import gzip
import argparse
import statistics as sts
from typing import List, Tuple

import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqFeature import SeqFeature
from Bio.SeqRecord import SeqRecord


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files
parser.add_argument(
    '-i',
    '--assm-acc-file',
    help="""TSV file (with header) with
  Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs""",
    required=True
)

parser.add_argument(
    '-g',
    '--gbk-dir',
    help='directory that contains downloaded gbk.gz files',
    required=True
)

# Output files
parser.add_argument(
    '-o',
    '--out-fasta',
    help='output fasta file for genes sequences',
    required=True
)

parser.add_argument(
    '-s',
    '--out-stats',
    help='output TSV file for per-replicon statistics',
    required=True
)

# Dependencies
parser.add_argument(
    '--cmsearch',
    help='cmsearch executable',
    required=True
)

parser.add_argument(
    '-r',
    '--rfam-family-cm',
    help=""".cm file containing covariance model of target gene family
  (RF00177 for bacterial ribosomal SSU, RF01959 for archaeal ribosomal SSU)""",
    required=True
)

parser.add_argument(
    '--seqkit',
    help='seqkit executable',
    required=True
)


args = parser.parse_args()


assm_acc_fpath = os.path.realpath(args.assm_acc_file)
gbk_dpath = os.path.realpath(args.gbk_dir)
fasta_outfpath = os.path.realpath(args.out_fasta)
outstats_fpath = os.path.realpath(args.out_stats)
cmsearch_fpath = os.path.realpath(args.cmsearch)
rfam_family_fpath = os.path.realpath(args.rfam_family_cm)
seqkit_fpath = os.path.realpath(args.seqkit)


# Check existance of input file -i/--assm-acc-file
if not os.path.exists(assm_acc_fpath):
    print(f'Error: file `{assm_acc_fpath}` does not exist!')
    sys.exit(1)
# end if

# Check existance of gbk directory -g/--gbk-dir
if not os.path.isdir(gbk_dpath):
    print(f'Error: directory `{gbk_dpath}` does not exist!')
    sys.exit(1)
# end if

# Create output directories if needed
for some_dir in map(os.path.dirname, [fasta_outfpath, outstats_fpath]):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end if

# Check existance of cmsearch executables --cmsearch and --seqkit
# And check if they are executable
for fpath in (cmsearch_fpath, seqkit_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
    # Check if an executable is actually executable
    if not os.access(fpath, os.X_OK):
        print(f'Error: file `{fpath}` is not executable!')
        sys.exit(1)
    # end if
# end for

# Check existance of file rfam_family_fpath --rfam-family-cm
if not os.path.exists(rfam_family_fpath):
    print(f'Error: file `{rfam_family_fpath}` does not exist!')
    sys.exit(1)
# end if


print(assm_acc_fpath)
print(gbk_dpath)
print(cmsearch_fpath)
print(rfam_family_fpath)
print(seqkit_fpath)
print()


# Header of cmsearch's .tblout output files
tblout_header = 'target_name\taccession\tquery_name\taccession\tmdl\tmdl_from\tmdl_to\tseq_from\tseq_to\tstrand\ttrunc\tpass\tgc\tbias\tscore\tEvalue\tinc\tdescription_of_target'

# Header of statistics file
stats_header = [
    'ass_id', 'gi_number', 'acc', 'title',
    'seq_start_truncation', 'improper_16S_annotation', 'topology',
    'num_genes', 'min_len', 'max_len', 'mean_len', 'median_len',
]

# Possible values of `product` qualifier if a 16S rRNA gene feature
ssu_product_names = {
    '16S RIBOSOMAL RNA',
    'SMALL SUBUNIT RIBOSOMAL RNA',
    'SMALL RIBOSOMAL RNA',
    'S-RRNA',
    'RIBOSOMAL RNA-16S'
}

# Notes which are assigned to truncated 16S genes
trunc_ssu_notes = {
    '16S RIBOSOMAL RNA RRNA PREDICTION IS TOO SHORT',
    'POSSIBLE 16S RIBOSOMAL RNA',
}


def is_ssu(feature: SeqFeature):
    # Script checks if feature is a 16S rRNA gene
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


def filter_ssu_genes(features: List[SeqFeature]):
    # Function filters features: it keeps only 16S rRNA genes
    return tuple(
        filter(
            is_ssu,
            features
        )
    )
# end def filter_ssu_genes


def is_annotated_with_pgap(gbrecord: SeqRecord):
    # Function checks if GenBank record is annotated with PGAP
    #   (NCBI Prokaryotic Genome Annotation Pipeline).
    # If it is -- returns True. If not -- returns False.

    # Find structured comment
    try:
        struct_comment = gbrecord.annotations['structured_comment']
    except KeyError:
        return False
    # end try

    # Find section "Genome-Annotation-Data" in structured comment
    if 'Genome-Annotation-Data' in struct_comment.keys():
        assembly_key = 'Genome-Annotation-Data'
    else:
        return False
    # end if

    # Get annotation data
    annotation_data = struct_comment[assembly_key]

    # Find field "Annotation Pipeline" in `annotation_data`
    if 'Annotation Pipeline' in annotation_data.keys():
        annot_pipe_key = 'Annotation Pipeline'
    else:
        return False
    # end if

    # Get field containing info about annotation pipeline
    annot_pipe = annotation_data[annot_pipe_key]

    # Find "NCBI PROKARYOTIC GENOME ANNOTATION PIPELINE" in this field
    if 'NCBI PROKARYOTIC GENOME ANNOTATION PIPELINE' in annot_pipe.upper():
        return True
    else:
        return False
    # end if
# end def


def seq_start_may_truncate_ssu(gbrecord: SeqRecord, features: List[SeqFeature]):
    # Function checks if gbrecord can contain 16S genes truncated by sequences start.
    # It it can -- returns True. therwise -- False.

    for f in features:
        if f.location.start == 0 or f.location.end == len(gbrecord.seq):
            return True
        # end if
    # end for
    return False
# end def seq_start_may_truncate_ssu


def extract_gene_as_is(feature: SeqFeature, gbrecord: SeqRecord):
    # Function extracts sequence of gene feature `feature` from `gbrecord`

    # We will report 1-based left-closed and right-closed coordinates
    seq_start = feature.location.start + 1
    seq_end = feature.location.end
    seq_strand = feature.location.strand

    # Extract sequence
    seq = gbrecord.seq[seq_start-1 : seq_end]

    # Make it reverse-complement if appropriate
    if seq_strand == -1:
        seq = seq.reverse_complement()
    # end if

    # Set `strand` value
    strand_str = 'plus' if seq_strand == 1 else 'minus'

    # Remove chars '<' and '>': they will appear if
    #   Bio.SeqFeature.BeforePosition or Bio.SeqFeature.AfterPosition
    #   are converted to str
    seq_start_for_header = str(seq_start).replace('<', '')
    seq_end_for_header = str(seq_end).replace('>', '')

    header = f'{gbrecord.id}:{seq_start_for_header}-{seq_end_for_header}_{strand_str} {gbrecord.description}'

    return header, seq
# end def extract_gene_as_is


def run_cmsearch(fasta_fpath: str):
    # Function runs cmsearch searching for 16S rRNA genes in sequence
    #   stored in file `fasta_fpath`
    # Returns path to result .tblout file

    tblout_fpath = '/tmp/tmpXXX_tblout.tsv'
    out_fpath = '/tmp/tmpXXX_cmsearch_out.txt'
    cmd = f'{cmsearch_fpath} --noali -o {out_fpath} --tblout {tblout_fpath} --cpu 6 {rfam_family_fpath} {fasta_fpath}'
    # cmd = f'{cmsearch_fpath} -o {out_fpath} --tblout {tblout_fpath} --cpu 6 {rfam_family_fpath} {fasta_fpath}'

    exit_code = os.system(cmd)
    if exit_code != 0:
        raise OSError(f'exit code {exit_code}')
    # end if
    return tblout_fpath
# end def run_cmsearch


def reformat_tblout(tblout_fpath: str):
    # Function makes raw file `tblout_fpath` (it's output of cmsearch) readable for pandas.read_csv

    # Read all lines not commented with #'s
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

    # Replace multiple spaces with single spaces
    for i in range(len(lines)):
        for space_num in range(20, 1, -1):
            lines[i] = lines[i].replace(' '*space_num, ' ')
        # end for
    # end for

    # Replace spaces in some vulnerable column values with underscores
    for i in range(len(lines)):
        lines[i] = lines[i].replace('Bacterial small subunit ribosomal RNA', 'Bacterial_small_subunit_ribosomal_RNA')
        lines[i] = lines[i].replace('Archaeal small subunit ribosomal RNA', 'Archaeal_small_subunit_ribosomal_RNA')
    # end for

    # Replace spaces with tabs
    for i in range(len(lines)):
        lines[i] = lines[i].replace(' ', '\t')
    # end for

    # Write result to the same .tblout file
    with open(tblout_fpath, 'wt') as tblout_file:
        tblout_file.write(f'{tblout_header}\n')
        tblout_file.write('\n'.join(lines) + '\n')
    # end with
# end def reformat_tblout


def amend_coordinates_on_extended_seq(tblout_df: pd.DataFrame, original_len: int) -> pd.DataFrame:
    # Function amends coordinates of SSU genes if there are any on extended tail

    def amend_length(row):
        # Amend coordinates of fixed copies of genes truncated by sequence start
        #   in order to replace coordinates greater than `original_len`
        if row['seq_from'] > original_len:
            row['seq_from'] -= original_len
        # end if
        if row['seq_to'] > original_len:
            row['seq_to'] -= original_len
        # end if
        return row
    # end def amend_length

    tblout_df = tblout_df.apply(amend_length, axis=1)
    return tblout_df
# end def amend_coordinates_on_extended_seq


def remove_sestart_truncated_gene(tblout_df: pd.DataFrame) -> pd.DataFrame:
    # Function removes from `tblout_df` genes truncated by sequence start if fixed variant
    #   of this gene exists, thanks to extended tail.

    def set_remove_flag(row):
        if row['strand'] == '+':
            if row['seq_from'] == 1: # if gene may be truncated by sequence start
                fixed_copy_exists = tblout_df[tblout_df['seq_to'] == row['seq_to']].shape[0] > 1
                if fixed_copy_exists: # if fixed copy of this gene exists in `tblout_df`
                    row['remove'] = 1
                # end if
            # end if
        else:
            if row['seq_to'] == 1: # if gene may be truncated by sequence start
                fixed_copy_exists = tblout_df[tblout_df['seq_from'] == row['seq_from']].shape[0] > 1
                if fixed_copy_exists: # if fixed copy of this gene exists in `tblout_df`
                    row['remove'] = 1
                # end if
            # end if
        # end if
        return row
    # end def set_remove_flag

    # Set flag `remove` to 1 on those genes, which should be removed
    tblout_df['remove'] = np.repeat(0, tblout_df.shape[0])
    tblout_df = tblout_df.apply(set_remove_flag, axis=1)

    # Remove truncated genes
    tblout_df = tblout_df[tblout_df['remove'] == 0]
    # tblout_df.drop(['remove'], axis=1, inplace=True)

    return tblout_df
# end def remove_sestart_truncated_gene


def extract_gene_seq_after_cmsearch(
    gbrecord: SeqRecord,
    seq_start: int,
    seq_end: int,
    topology: str) -> Seq:
    # Function extracts gene sequence from `gbrecord`.
    # `seq_start` and `seq_end` are 1-based, left-closed, right-closed.
    # Function performs well even if `seq_start` > `seq_end`: it may occur if
    #   gene is truncated by sequence start.

    if seq_start < seq_end:
        seq = gbrecord.seq[seq_start-1 : seq_end]
    elif seq_start > seq_end and topology == 'circular':
        seq = gbrecord.seq[seq_start-1 :] + gbrecord.seq[: seq_end]
    else:
        print('\nError 2!')
        print('After cmsearch, seq_start > seq_end and topology is not "circular"')
        print('Please, contact the developer: he thought, this case is impossible.')
        print(f"""Debug info:
    ACC: {gbrecord.id};
    seq_start={seq_start};
    seq_end={seq_end};
    topology=`{topology}`""")
        sys.exit(2)
    # end if

    return seq
# end def extract_gene_seq_1based_coords


def extract_reannotated_genes(gbrecord: SeqRecord, topology: str):
    # Function reannotates 16S rRNA genes in `gbrecord` with cmsearch
    #   and extracts sequences of discovered genes from it.

    # Save original length of sequence
    original_len = len(gbrecord.seq)
    len_circ_tail = min(original_len, 2000)

    # Make circular sequence a bit longer in order to annotate properly
    #   genes truncated by sequecne start
    if topology == 'circular':
        gbrecord.seq = gbrecord.seq + gbrecord.seq[: len_circ_tail]
    # end if

    # Perform reannotations
    tmp_fasta = f'tmp{os.getpid()}.fasta'
    with open(tmp_fasta, 'w') as tmpf:
        tmpf.write(f'>{gbrecord.id}\n{str(gbrecord.seq)}\n')
    # end with
    tblout_fpath = run_cmsearch(tmp_fasta)
    reformat_tblout(tblout_fpath)
    os.unlink(tmp_fasta)

    # Now we have `tblout_fpath` and cam extract genes sequences from `gbrecord`

    tblout_df = pd.read_csv(tblout_fpath, sep='\t')

    if topology == 'circular':
        tblout_df = amend_coordinates_on_extended_seq(
            tblout_df,
            original_len
        )
        tblout_df = remove_sestart_truncated_gene(tblout_df)
        gbrecord.seq = gbrecord.seq[: -len_circ_tail]
    # end if


    genes = list()

    # Iterate over rows of tblout_df and extract sequences of annotated genes
    for _, row in tblout_df.iterrows():
        seq_start = row['seq_from']
        seq_end = row['seq_to']
        seq_strand = row['strand']

        if seq_strand == '+':
            seq = extract_gene_seq_after_cmsearch(gbrecord, seq_start, seq_end, topology)
        else:
            seq_start, seq_end = seq_end, seq_start
            seq = extract_gene_seq_after_cmsearch(gbrecord, seq_start, seq_end, topology)
            seq = seq.reverse_complement()
        # end if

        strand_str = 'plus' if seq_strand == '+' else 'minus'

        seq_start_for_header = seq_start
        seq_end_for_header = seq_end

        # Configure sequence header
        seq_header = '{}:{}-{}_{} {}'\
            .format(
                gbrecord.id,
                seq_start_for_header,
                seq_end_for_header,
                strand_str,
                gbrecord.description
            )

        # Add extracted gene to the list
        genes.append(
            (
                seq_header,
                seq
            )
        )
    # end for

    return genes
# end def extract_reannotated_genes


def get_seq_lengths(extracted_genes: List[Tuple[str, str]]):
    # Function returns lengths of genes passed to it.
    # `extracted_genes` is an object retuned by `extract_reannotated_genes`
    lengths = [len(gene[1]) for gene in extracted_genes]
    return lengths
# end def get_seq_lengths


def calc_gene_stats(extracted_genes: List[Tuple[str, str]]):
    # Function calculates statistics for gene set passed to it:
    #  number of genes
    #  min length
    #  max length
    #  mean length
    #  median length
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


# Read input data

acc_df = pd.read_csv(
    assm_acc_fpath,
    sep='\t'
)


# accs_select = {
#     'NZ_CP013210.1',
#     'NZ_CP008696.1',
#     # 'NZ_CP050525.1',
# }
# acc_df = acc_df.query('acc in @accs_select')

n_accs = acc_df.shape[0] # number of ACCESSION.VESRION's to process



# == Proceed ==

with open(fasta_outfpath, 'wt') as fasta_outfile, open(outstats_fpath, 'wt') as stats_outfile:

    # Write header of stats file
    stats_outfile.write('{}\n'.format('\t'.join(stats_header)))

    # For each RefSeq record: extract 16S rRNA genes from it
    for i, row in acc_df.iterrows():

        ass_id = row['ass_id']
        gi_number = row['gi_number']
        acc = row['acc']
        title = row['title']

        print(f'\rDoing {i+1}/{n_accs}: {acc}', end=' '*10)

        # Configure output file path (e.g. `NZ_CP063178.1.gbk.gz`)
        gbk_fpath = os.path.join(
            gbk_dpath,
            f'{acc}.gbk.gz'
        )

        # Read GenBank file
        with gzip.open(gbk_fpath, 'rt') as gbfile:
            gbrecord = tuple(SeqIO.parse(gbfile, 'gb'))[0]
        # end with

        # Get 16S rRNA features
        ssu_features = filter_ssu_genes(gbrecord.features)

        # Check if 16S rRNA genes may be truncated by sequence start
        seq_start_truncation = seq_start_may_truncate_ssu(gbrecord, ssu_features)

        # Check if sequenec was annotated with PGAP
        improper_16S_annotation = not is_annotated_with_pgap(gbrecord)

        # Get topology
        topology = None
        if 'topology' in gbrecord.annotations.keys():
            topology = gbrecord.annotations['topology']
        # end if

        if improper_16S_annotation or (seq_start_truncation and topology == 'circular'):
            # Reannotate 16S rRNA genes using cmsearch
            # And extract reannotated genes
            extracted_genes = extract_reannotated_genes(gbrecord, topology)
            for header, seq in extracted_genes:
                fasta_outfile.write(f'>{header}\n{seq}\n')
            # end for
        else:
            # Just extract genes, that have been already annotated
            extracted_genes = [extract_gene_as_is(f, gbrecord) for f in ssu_features]
            for header, seq in extracted_genes:
                fasta_outfile.write(f'>{header}\n{seq}\n')
            # end for
        # end if

        # Write statistics to stats file
        stats_outfile.write(f'{ass_id}\t{gi_number}\t{acc}\t{title}\t')
        stats_outfile.write(f'{1 if seq_start_truncation else 0}\t{1 if improper_16S_annotation else 0}\t{topology}\t')

        num_genes, min_len, max_len, mean_len, median_len = calc_gene_stats(extracted_genes)
        stats_outfile.write(f'{num_genes}\t{min_len}\t{max_len}\t{mean_len}\t{median_len}\n')
    # end for
# end with


# Some genes may be duplicated because of appending sequence's start to it's end
# Therefore, we need to dereplicate sequences by name (seqkit rmdup -n)
print('\nRunning seqkit rmdup...')
tmpfasta = os.path.join(
    os.path.dirname(fasta_outfpath),
    'tmp_ALL_GENES.fasta'
)
os.system(f'cat {fasta_outfpath} | {seqkit_fpath} rmdup -n > {tmpfasta}')
os.system(f'cat {tmpfasta} | {seqkit_fpath} seq -u > {fasta_outfpath}')
os.unlink(tmpfasta)
print('Done\n')

print('\nCompleted!')
print(fasta_outfpath)
print(outstats_fpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
