#!/usr/bin/env python3

# TODO: update description

# The script assigns categories to downlaoded genomes. Categories are assigned according
#   to the reliability of a genome assembly.
# The categories are the following:
# Category 1. A genome is not of category 3, and it was sequenced using PacBio or ONT+Illumina.
# Category 2. A genome is not of category 3, and it was sequenced neither using PacBio nor ONT+Illumina.
# Category 3. At least one of the following is true:
#   - A genome has at least one degenerate base in its SSU gene sequences.
#   - At least one of the genomic sequences contains phrase "map unlocalized" in it title,
#     and the sequence contains an SSU gene (or a part of it).

## Command line arguments

### Input files:
# 1. `-f / --all-fasta-file` -- a fasta file with all extracted genes sequences.
#   This file is the output of the script `extract_16S.py`. Mandatory.
# 2. `-s / --all-stats-file` -- a file with per-replicon statistics of input 16S gene sequences.
#   This file is the output of the script `extract_16S.py`, too. Mandatory.
# 3. `-g / --gbk-dir` -- the directory where the downloaded `.gbk.gz` files are located
#   (see script `download_genomes.py`). Mandatory.

### Output files:
# 1. `-o / --outfile` -- an output file mapping seqIDs (and Assembly IDs) to categories. Mandatory.
# 2. `-l / --seqtech-logfile` -- a log file to track if sequence technology is successfully
#   extracted from all genomes where it is specified. "Seqtech" means SEQuencing TECHnology. Mandatory.

### Dependencies:
# 1. `--seqkit` -- a `seqkit` executable: github.com/shenwei356/seqkit. Mandatory.


import os

print(f'\n|=== STARTING SCRIPT `{os.path.basename(__file__)}` ===|\n')


import sys
import gzip
import argparse
import subprocess as sp
from typing import Dict, List, Sequence, TextIO

import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--all-fasta-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-s',
    '--all-stats-file',
    help="""TSV file (with header) containing per-replicons SSU gene statistics
    reported by extract_16S.py""",
    required=True
)

parser.add_argument(
    '-g',
    '--gbk-dir',
    help='directory that contains downloaded gbk.gz files',
    required=True
)

parser.add_argument(
    '--prev-categories',
    help='file of genome categories from previous RiboGrove release',
    required=False
)
parser.add_argument(
    '--prev-assm-acc-file',
    help="""TSV file (with header) with
    Assembly IDs, GI numbers, ACCESSION.VERSION's and titles separated by tabs
    from previous RiboGrove release""",
    required=False
)


# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='output file mapping genes seqIDs and Assembly IDs to categories',
    required=True
)

parser.add_argument(
    '-l',
    '--seqtech-logfile',
    help="""log file to track if sequencing technology is successfully extracted
    from all genomes where it is specified""",
    required=True
)

# Dependencies

parser.add_argument(
    '--seqkit',
    help='seqkit executable',
    required=True
)


args = parser.parse_args()

# For convenience
fasta_seqs_fpath = os.path.abspath(args.all_fasta_file)
in_stats_fpath = os.path.abspath(args.all_stats_file)
gbk_dpath = os.path.abspath(args.gbk_dir)
if not args.prev_categories is None and not args.prev_assm_acc_file is None:
    cached_categories = True
    prev_categories_fpath = os.path.abspath(args.prev_categories)
    prev_assm_acc_fpath = os.path.abspath(args.prev_assm_acc_file)
else:
    cached_categories = False
    prev_categories_fpath = None
    prev_assm_acc_fpath = None
# end if
outfpath = os.path.abspath(args.outfile)
seqtech_logfpath = os.path.abspath(args.seqtech_logfile)
seqkit_fpath = os.path.abspath(args.seqkit)


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, in_stats_fpath, seqkit_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Check if input directory exists
if not os.path.exists(gbk_dpath):
    print(f'Error: directory `{gbk_dpath}` does not exist!')
    sys.exit(1)
# end if

# Check if seqkit executable is actually executable
if not os.access(seqkit_fpath, os.X_OK):
    print(f'Error: file `{seqkit_fpath}` is not executable!')
    sys.exit(1)
# end if

# Create output directories if needed
for some_dir in map(os.path.dirname, [outfpath, seqtech_logfpath]):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end if

# Check if previous ("cached") files is specified
if cached_categories:
    for f in (prev_categories_fpath, prev_assm_acc_fpath):
        if not os.path.exists(f):
            print(f'Error: file `{f}` does not exist!')
            sys.exit(1)
        # end if
    # end for
# end if

print(fasta_seqs_fpath)
print(in_stats_fpath)
print(gbk_dpath)
print(prev_categories_fpath)
print(prev_assm_acc_fpath)
print(seqtech_logfpath)
print(seqkit_fpath)
print()


# Paths to files containing merker words fow identifying sequencing technologies
pacbio_vocab_fpath = os.path.realpath(
    os.path.join(os.path.dirname(__file__), 'seqtech_dicts/pacbio')
)
illumina_vocab_fpath = os.path.realpath(
    os.path.join(os.path.dirname(__file__), 'seqtech_dicts/illumina')
)
nanopore_vocab_fpath = os.path.realpath(
    os.path.join(os.path.dirname(__file__), 'seqtech_dicts/ont')
)


# These keywords might be recognized as "ONT", but they are not
#   releated to Oxford Nanopore.
# These words must be removed from seqtech string befor searching for keywords
seem_like_ont_but_not = {
    'IONTORRENT',
    # ~~~
    'CONTIG', # see NC_020549.1
    # ~~~
}



def find_degenerate_in_16S(fasta_seqs_fpath: str, stats_df: pd.DataFrame, seqkit_fpath: str) -> List[str]:
    # Function reports Assembly IDs of genomes which contain degenerate bases in their SSU genes.

    # Configure command reporting accessions of sequences which contain degenerate bases in their SSU genes
    cmd = f'{seqkit_fpath} grep -srp "[RYWSKMHVBDN]" {fasta_seqs_fpath} | {seqkit_fpath} seq -ni | cut -f1 -d":" | sort | uniq'
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error at popen finding degenerate in 16S genes')
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(1)
    else:
        # Parse accessions
        accs_degenerate_in_16S = set(stdout_stderr[0].decode('utf-8').split('\n'))
    # end if

    # Translate obtained accessions to Assembly IDs
    ass_ids_degenerate_in_16S = set(stats_df.query('acc in @accs_degenerate_in_16S')['ass_id'])

    return ass_ids_degenerate_in_16S
# end def find_degenerate_in_16S


# TODO: remove seqkit from here
def get_genes_seqIDs(fasta_seqs_fpath: str, seqkit_fpath: str) -> List[str]:
    # Function reports all seqIDs of sequences from given fasta file fasta_seqs_fpath.

    # Configure command reporting seqIDs of fasta file
    cmd = f'{seqkit_fpath} seq -ni {fasta_seqs_fpath}'
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        print('Error at popen: extracting genes\' seqIDs')
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(pipe.returncode)
    else:
        # Parse seqIDs
        genes_seqIDs = list(stdout_stderr[0].decode('utf-8').split('\n'))
    # end if

    return genes_seqIDs
# end def get_genes_seqIDs


def make_acc_seqIDs_dict(fasta_seqs_fpath: str, seqkit_fpath: str) -> Dict[str, List[str]]:
    # Function creates dictionary that maps accessions to seqIDs

    # Get all seqIDs of gene sequences.
    # We want to have our seqIDs in the same order as they exist in original fasta file,
    #   therefore, we do this `reversed` here
    genes_seqIDs = list(
        reversed(
            get_genes_seqIDs(fasta_seqs_fpath, seqkit_fpath)
        )
    )

    acc_seqIDs_dict = dict()

    for _ in range(len(genes_seqIDs)):

        seqID = genes_seqIDs.pop() # get next seqID
        acc = seqID.partition(':')[0] # parse ACCESSION.VERSION from seqID

        # Fill dictionary
        try:
            acc_seqIDs_dict[acc].append(seqID)
        except KeyError:
            acc_seqIDs_dict[acc] = [seqID]
        # end try
    # end for

    return acc_seqIDs_dict
# end def make_acc_seqIDs_dict


def read_seqtech_vocab(fpath: str) -> Sequence[str]:
    # Function makes vocabulary of marker seqtech words
    with open(fpath, 'rt') as vocab_file:
        vocab = tuple(
            map(
                lambda s: s.strip().upper(),
                vocab_file.readlines()
            )
        )
    # end with
    return vocab
# end


def is_pacbio(seqtech_str: str) -> bool:
    # Function decides if genome was sequenced with PacBio

    global pacbio_vocab

    return any(
        tuple(
            map(
                lambda x: x in seqtech_str,
                pacbio_vocab
            )
        )
    )
# end

def is_illumina(seqtech_str: str) -> bool:
    # Function decides if genome was sequenced with Illumina

    global illumina_vocab

    return any(
        tuple(
            map(
                lambda x: x in seqtech_str,
                illumina_vocab
            )
        )
    )
# end

def is_nanopore(seqtech_str: str) -> bool:
    # Function decides if genome was sequenced with Oxford Nanopore

    global nanopore_vocab

    for keyword in nanopore_vocab:
        if keyword != 'ONT':
            # If keyword is not "ONT" we proceed just like with any other keyword
            if keyword in seqtech_str:
                return True
            # end if
        else:
            # It kwyword is "ONT", we should previously remove words from `seem_like_ont_but_not`
            #   before searching the keyword
            for word in seem_like_ont_but_not:
                seqtech_str = seqtech_str.replace(word, '')
            # end for
            if keyword in seqtech_str:
                return True
            # end if
        # end if
    # end for
    return False
# end def is_nanopore


def parse_seqtech(gbrecord: SeqRecord, logfile: TextIO) -> str:
    # Function parses Sequencing Technology from sequence record (gbrecord)

    # Find structural comment
    try:
        struct_comment = gbrecord.annotations['structured_comment']
    except KeyError as err:
        logfile.write(f'{gbrecord.id} - Error (no structured_comment): {err}\n')
        return None
    # end try

    # Find fiels "Genome-Assembly-Data" or "Assembly-Data" in the structured comment
    if 'Genome-Assembly-Data' in struct_comment.keys():
        assembly_key = 'Genome-Assembly-Data'
    elif 'Assembly-Data' in struct_comment.keys():
        assembly_key = 'Assembly-Data'
    else:
        logfile.write(f'{gbrecord.id} - Error: no `(Genome)-Assembly-Data` in keys of structured_comment. ')
        logfile.write(f'Keys: {";".join(struct_comment.keys())}\n')
        return None
    # end if

    # Get assembly data
    assembly_data = struct_comment[assembly_key]

    # Find key for extracting sequencing technology
    if 'Sequencing Technology' in assembly_data.keys():
        seqtech_key = 'Sequencing Technology'
    elif 'Sequencing Technolog' in assembly_data.keys():
        seqtech_key = 'Sequencing Technolog'
    elif 'Sequencing technology' in assembly_data.keys():
        seqtech_key = 'Sequencing technology'
    else:
        logfile.write(f'{gbrecord.id} - Error: no `Sequencing Technology` in keys of `Assembly data` ')
        logfile.write(f'Keys: {";".join(assembly_data.keys())}\n')
        return None
    # end if

    # Get sequencing technology
    seqtech = assembly_data[seqtech_key]

    # Save result in the log file
    logfile.write(f'{gbrecord.id} - ok\n')

    return seqtech.upper().strip()
# end def parse_seqtech


def find_NNN(gbrecord: SeqRecord) -> bool:
    # Function checks if there is an 'NNN' substring in sequence of `gbrecord`
    return 'NNN' in str(gbrecord.seq)
# end def find_NNN


# Read per-replicon statistics
stats_df = pd.read_csv(
    in_stats_fpath,
    sep='\t'
)
n_accs = stats_df.shape[0] # count replicon sequences


# Make vocabularies of seqtech keywords
pacbio_vocab = read_seqtech_vocab(pacbio_vocab_fpath)
illumina_vocab = read_seqtech_vocab(illumina_vocab_fpath)
nanopore_vocab = read_seqtech_vocab(nanopore_vocab_fpath)


# Get Assembly IDs of genomes having degenerate bases in SSU genes
print('Searching for 16S genes containing degenerate bases...')
ass_ids_degenerate_in_16S = find_degenerate_in_16S(fasta_seqs_fpath, stats_df, seqkit_fpath)
print(f'Found {len(ass_ids_degenerate_in_16S)} assemblies containing 16S genes with degenerate bases')

# Create dictionary that maps ACCESSION.VERSION's to seqIDs
print('Building auxiliary data structures...')
acc_seqIDs_dict = make_acc_seqIDs_dict(fasta_seqs_fpath, seqkit_fpath)
print('done')


def encode_accs(acc_list):
    return ''.join(sorted(acc_list))
# end def

if cached_categories:
    print('Reading cached files')

    prev_categories_df = pd.read_csv(prev_categories_fpath, sep='\t')
    cached_ass_ids = set(prev_categories_df['ass_id'])

    prev_assm_acc_df = pd.read_csv(prev_assm_acc_fpath, sep='\t')
    ass_ids = set(prev_assm_acc_df['ass_id'])
    acc_code_dict = dict()
    for ass_id in ass_ids:
        accs = tuple(prev_assm_acc_df[prev_assm_acc_df['ass_id'] == ass_id]['acc'])
        acc_code_dict[ass_id] = encode_accs(accs)
    # end for
    print('done')
else:
    cached_ass_ids = set()
# end if


# == Proceed ==

with open(outfpath, 'wt') as outfile, \
     open(seqtech_logfpath, 'wt') as logfile:

     # Write headers
    outfile.write('ass_id\tseqID\tseqtech\tcontains_NNN\tdegenerate_in_16S\tunlocalized_16S\tcategory\n')

    # Get all Assembly IDs
    assembly_IDs = tuple(set(stats_df['ass_id']))

    # Iterate over Assembly IDs
    for i, ass_id in enumerate(assembly_IDs):
        print(f'\rDoing {i+1}/{len(assembly_IDs)}: {ass_id}', end=' '*10)

        # Try find cached data
        if ass_id in cached_ass_ids:
            accs = tuple(stats_df[stats_df['ass_id'] == ass_id]['acc'])
            acc_code = encode_accs(accs)
            if acc_code == acc_code_dict[ass_id]:
                cached_ass_df = prev_categories_df[prev_categories_df['ass_id'] == ass_id]
                cached_ass_df.to_csv(
                    outfile,
                    sep='\t',
                    header=False,
                    index=False,
                    na_rep='NA',
                    encoding='utf-8'
                )
                continue # cache hit
            # end if
        # end if

        # Get rows corresponding to current assembly
        ass_df = stats_df[stats_df['ass_id'] == ass_id]


        # Genome has (maybe, patrial) SSU genes in "map unlocalized" sequences
        unlocalized_16S = False

        # Genome contains at least 3 N's in a row
        contains_NNN = False

        # Genome has degenerate bases in SSU genes
        degenerate_in_16S = ass_id in ass_ids_degenerate_in_16S

        # List of seqtech strings
        seqtechs = list()

        # Iterate over ACCESSION.VERSION's of current genome
        for acc_i, row in ass_df.iterrows():
            acc = row['acc']
            title = row['title']

            # Update flag `unlocalized_16S`
            map_unlocalized = 'MAP UNLOCALIZED' in title.upper()
            unlocalized_16S = unlocalized_16S or (map_unlocalized and row['num_genes'] != 0)

            # Configure path to GenBank file containing current sequence
            gbk_fpath = os.path.join(
                gbk_dpath,
                f'{acc}.gbk.gz'
            )

            # Read GenBank record
            with gzip.open(gbk_fpath, 'rt') as gbfile:
                gbrecord = tuple(SeqIO.parse(gbfile, 'gb'))[0]
            # end with

            # Update flag `contains_NNN`
            contains_NNN = contains_NNN or find_NNN(gbrecord)

            # Parse setech and add it to `seqtechs`
            curr_seqtech = parse_seqtech(gbrecord, logfile)
            if not curr_seqtech is None:
                seqtechs.append(curr_seqtech)
            # end if
        # end for

        if len(seqtechs) != 0:
            # Seqtech is found and we will merely concatenate them
            seqtech = '. '.join(seqtechs)
        else:
            seqtech = 'NA'
        # end if

        # Assign category to the genome
        category = None
        if contains_NNN:
            category = 'NA'
        elif degenerate_in_16S or unlocalized_16S:
            category = 3
        elif is_pacbio(seqtech) or (is_illumina(seqtech) and is_nanopore(seqtech)):
            category = 1
        else:
            category = 2
        # end if

        # Get all ACCESSION.VERSION's of current genome
        accs = tuple(ass_df['acc'])

        # Write ouput to output file
        for acc in accs:
            try:
                for seqID in acc_seqIDs_dict[acc]:
                    outfile.write(f'{ass_id}\t{seqID}\t{seqtech}\t')
                    outfile.write(f'{1 if contains_NNN else 0}\t{1 if degenerate_in_16S else 0}\t')
                    outfile.write(f'{1 if unlocalized_16S else 0}\t{category}\n')
                # end for
            except KeyError:
                pass
            # end try
        # end for
    # end for
# end with

print('\nCompleted!')
print(outfpath)
print(seqtech_logfpath)
print(f'\n|=== EXITTING SCRIPT `{os.path.basename(__file__)}` ===|\n')
