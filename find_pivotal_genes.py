#!/usr/bin/env python3

# Script finds so called pivotal genes. Pivotal gene (within a genome) is SSU gene,
#   which is most likely functional, i.e. it has the highest score when
#   aligning Rfam SSU covariance model against all SSU genes sequenes from this genome.

# Input files:
# 1. Fasta file of genes sequences (-f/--fasta-seqs-file).
# 2. TSV file if per-replicon genes statistics (-s/--genes-stats-file).

# Output files:
# 1. TSV file reporting pivotal genes for each genome, in which maximum length difference
#    between SSU genes is greater than `lendiff_threshold`.
# 2. Directory for storing .tblout files -- output files pf cmscan (-t/--tblout-dir).

# Dependencies:
# 1. cmscan from Infernal: http://eddylab.org/infernal/ (--cmscan)
# 2. cmpress from Infernal: http://eddylab.org/infernal/ (--cmpress)
# 3. .cm file containing covariance model of target gene family (--rfam-family-cm)
#   (RF00177 for bacterial ribosomal SSU, RF01959 for archaeal ribosomal SSU)

# Parameters:
# 1. --lendiff-threshold
#    It is the length difference between SSU genes in a genome to start searching for pivotal gene(s).
#    If maxinum length difference between SSU genes in a genome is greter than `--lendiff-threshold`,
#    the script will search for pivotal genes, i.e. it will launch cmscan.


import os
import re
import sys
import argparse
import subprocess as sp
from typing import Sequence, Tuple, List, Dict

import pandas as pd
from Bio import SeqIO



# == Parse arguments ==

parser = argparse.ArgumentParser()

# Input files

parser.add_argument(
    '-f',
    '--fasta-seqs-file',
    help='fasta file of SSU gene sequences',
    required=True
)

parser.add_argument(
    '-s',
    '--genes-stats-file',
    help='TSV file (with header) containing per-replicons SSU gene statistics',
    required=True
)

# Output files

parser.add_argument(
    '-o',
    '--outfile',
    help='output file reporting pivotal genes for each genome',
    required=True
)

parser.add_argument(
    '-t',
    '--tblout-dir',
    help='output directoryfor storing .tblout files -- output files of cmscan',
    required=True
)


# Dependencies

parser.add_argument(
    '--cmscan',
    help='cmscan executable',
    required=True
)

parser.add_argument(
    '--cmpress',
    help='cmpress executable',
    required=True
)

parser.add_argument(
    '--rfam-family-cm',
    help=""".cm file containing covariance model of target gene family
  (RF00177 for bacterial ribosomal SSU, RF01959 for archaeal ribosomal SSU)""",
    required=True
)

# Heuristic's params
parser.add_argument(
    '--lendiff-threshold',
    help="""Length difference between SSU genes in a genome to start searching for pivotal gene(s).
    If maxinum length difference between SSU genes in a genome is greter than `--lendiff-threshold`,
    the script will search for pivotal genes, i.e. it will launch cmscan.""",
    required=True
)


args = parser.parse_args()


# For convenience
fasta_seqs_fpath = os.path.abspath(args.fasta_seqs_file)
genes_stats_fpath = os.path.abspath(args.genes_stats_file)
outfpath = os.path.abspath(args.outfile)
tblout_dpath = os.path.abspath(args.tblout_dir)
cmscan_fpath = os.path.abspath(args.cmscan)
cmpress_fpath = os.path.abspath(args.cmpress)
rfam_fpath = os.path.abspath(args.rfam_family_cm)


# Check (and assign) value of `--lendiff-threshold`
try:
    lendiff_threshold = int(args.lendiff_threshold)
    if lendiff_threshold < 0:
        raise ValueError
    # end if
except ValueError:
    print('Error: you entered and invalid `--lendiff-threshold` value.')
    print(f'Your value: `{args.lendiff_threshold}`')
    print('It must be integer number >= 0')
    sys.exit(1)
# end try


# Check existance of all input files and dependencies
for fpath in (fasta_seqs_fpath, genes_stats_fpath, cmscan_fpath, cmpress_fpath, rfam_fpath):
    if not os.path.exists(fpath):
        print(f'Error: file `{fpath}` does not exist!')
        sys.exit(1)
    # end if
# enb for

# Check if executables are actually executable
for exetutable in (cmscan_fpath, cmpress_fpath):
    if not os.access(exetutable, os.X_OK):
        print(f'Error: file `{exetutable}` is not executable!')
        sys.exit(1)
    # end if
# end for

# Create output directory if needed
for some_dir in (tblout_dpath, os.path.dirname(outfpath)):
    if not os.path.isdir(some_dir):
        try:
            os.makedirs(some_dir)
        except OSError as err:
            print(f'Error: cannot create directory `{some_dir}`')
            sys.exit(1)
        # end try
    # end if
# end for


# Heder for reformatted .tblout files
tblout_header = 'target_name\taccession\tquery_name\taccession\tmdl\tmdl_from\tmdl_to\tseq_from\tseq_to\tstrand\ttrunc\tpass\tgc\tbias\tscore\tEvalue\tinc\tdescription_of_target'

# Temporary fasta file for storing queries for cmscan
query_fasta_fpath = 'tmpQUERY.fasta'


def run_cmpress(cmpress_fpath: str, rfam_fpath: str) -> None:
    # Function runs cmpress on rfam .cm file.
    # It is required to run cmscan.

    # Actually run cmpress
    print('Running `cmpress`...')
    cmd = f'{cmpress_fpath} {rfam_fpath}'
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        # It something goes wrong -- check the error message
        error_msg = stdout_stderr[1].decode('utf-8')

        # If `error_msg` contains `already_exists_msg_pattern` -- it's ok --
        #   index files exist.
        already_exists_msg_pattern = r'.+ file ({}.+) already exists'.format(rfam_fpath)
        already_exists_obj = re.search(already_exists_msg_pattern, error_msg)
        just_already_exists = not already_exists_obj is None

        if just_already_exists:
            print(error_msg)
            print(f'Removing {already_exists_obj.group(1)}')
            os.unlink(already_exists_obj.group(1))
            run_cmpress(cmpress_fpath, rfam_fpath)
        else:
            # If `error_msg` does not contain `already_exists_msg_pattern` -- oh, we must terminate
            print('Error: cannot cmpress .cm file')
            print(error_msg)
            sys.exit(1)
        # end if
    else:
        # Print piped stdout
        print(stdout_stderr[0].decode('utf-8'))
    # end if
# end def run_cmpress


def select_seqs(accs: Sequence, fasta_seqs_fpath: str, query_fasta_fpath: str) -> None:
    # Funciton selects sequences from file `fasta_seqs_fpath` and writes
    #   them to file `query_fasta_fpath`.
    # Function selects sequences assotiated with ACCESSION.VERSION's in tuple `accs`.

    # Configure command for selecting seuences from `accs` records
    acc_options = '-p "' + '" -p "'.join(accs) + '"'
    cmd = f'cat {fasta_seqs_fpath} | seqkit grep -nr {acc_options} > {query_fasta_fpath}'

    # Run command
    returncode = os.system(cmd)
    if returncode != 0:
        print(f'Error: cannot select sequences from `{fasta_seqs_fpath}`!')
        print(f'Accessions: {", ".join(accs)}')
        sys.exit(1)
    # end if
# end def select_seqs


def run_cmscan(
    cmscan_fpath: str,
    query_fasta_fpath: str,
    rfam_fpath: str,
    tblout_fpath: str,
    tblout_header) -> Tuple[str, pd.DataFrame]:
    # Function runs cmscan: it searches genes sequencs against a covariance model database (`rfam_fpath`).


    # Configure cpmmand
    cmd = f'{cmscan_fpath} --tblout {tblout_fpath} --toponly --cpu 6 --acc {rfam_fpath} {query_fasta_fpath}'

    # Run command
    pipe = sp.Popen(cmd, shell=True, stdout=sp.PIPE, stderr=sp.PIPE)
    stdout_stderr = pipe.communicate()

    if pipe.returncode != 0:
        # handle possible error
        print('Error at running cmscan!')
        print(stdout_stderr[1].decode('utf-8'))
        sys.exit(1)
    # end if

    # Parse output
    out_text = stdout_stderr[0].decode('utf-8').splitlines()

    # Reformat .tblout file: it has awful structure
    reformat_tblout(tblout_fpath, tblout_header)

    # Read DataFrame from reformatted .tblout file
    tblout_df = pd.read_csv(tblout_fpath, sep='\t')
    tblout_df.index = tblout_df['query_name'] # for convenience

    return out_text, tblout_df
# end def run_cmscan


def reformat_tblout(tblout_fpath: str, tblout_header: str) -> None:
    # Function reformats .tblout file, which is output of cmscan.
    # Raw .tblout file is awfully formatted: fields are delimited by spaces
    #   (sometimes 2, sometimes 3 etc). And columns names caontain spaces too.
    # We will reformat tis file: rename columns according to `tblout_header`
    #   and replace irregular field delimiter with tabs.

    # Read all lines except of those starting with #
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

    # Remove multiple spaces with single space
    for i in range(len(lines)):
        for space_num in range(20, 1, -1):
            lines[i] = lines[i].replace(' '*space_num, ' ')
        # end for
    # end for

    # Remove spaces in important fields with underscores
    for i in range(len(lines)):
        lines[i] = lines[i].replace('Bacterial small subunit ribosomal RNA', 'Bacterial_small_subunit_ribosomal_RNA')
        lines[i] = lines[i].replace('Archaeal small subunit ribosomal RNA', 'Archaeal_small_subunit_ribosomal_RNA')
    # end for

    # Replace spaces with tabs
    for i in range(len(lines)):
        lines[i] = lines[i].replace(' ', '\t')
    # end for

    # Write result lines to original file
    with open(tblout_fpath, 'wt') as tblout_file:
        tblout_file.write(f'{tblout_header}\n')
        tblout_file.write('\n'.join(lines) + '\n')
    # end with
# end def reformat_tblout


def parse_seqIDs(query_fasta_fpath: str) -> List[str]:
    # Function parses seqIDs from a fasta file.
    seq_records = SeqIO.parse(query_fasta_fpath, 'fasta')
    return list(map(lambda x: x.id, seq_records))
# end def parse_seqIDs


def amend_scores(seqIDs: List[str], out_text: str, tblout_df: pd.DataFrame) -> None:
    # Function amends scores calculated by cmscan.
    # The reason is that cmscan does not penalize some insertions (they call them "local end alignments").
    # But we do want to penalized these insertions, and this function does it.
    # The function modifies `tblout_df` -- this is the output.

    # Insert pattern: like *[183]*
    insrt_pattern = r'\*\[([0-9]+)\]\*'

    for seqID in seqIDs:

        if seqID in tblout_df.index:
            # Parse alienment: it contains out seqID in every alignment line
            curr_out_text = ''.join(filter(lambda l: seqID in l, out_text))

            # Find local end alignments
            unpenalted_insertions = re.findall(insrt_pattern, curr_out_text)

            # Calculate penalty
            penalty = sum(map(int, unpenalted_insertions))

            # Subtract penalty from initial score
            tblout_df.loc[seqID, 'score'] = tblout_df.loc[seqID, 'score'] - penalty
        # end if
    # end for
# end def amend_scores


def extract_pivotal_gene_lengths(query_fasta_fpath: str, best_gene_seqIDs: Sequence[str]) -> Dict[str, int]:
    # Function creates dictionary, which maps seqIDs of pivotal genes
    #   to lengths of corresponding sequences.

    # Read seq records
    seq_records = tuple(SeqIO.parse(query_fasta_fpath, 'fasta'))

    # Filter picotal genes
    pivotal_records = filter(
        lambda rec: rec.id in best_gene_seqIDs,
        seq_records
    )

    # Cofigure dictionary and return it
    return {r.id: len(r.seq) for r in pivotal_records}
# end def extract_pivotal_genes


# == Proceed ==

# Read statistics
stats_df = pd.read_csv(genes_stats_fpath, sep='\t')

# Transform per-replicon statistics to per-genome statistics
grpd_df = stats_df.groupby('ass_id').agg({'min_len': 'min', 'max_len': 'max'}).reset_index()
grpd_df['lendiff'] = grpd_df['max_len'] - grpd_df['min_len']

# Index file with covariance model
run_cmpress(cmpress_fpath, rfam_fpath)


print('Searching for pivotal genes...')

# Assembly IDs
ass_ids = set(stats_df['ass_id'])


with open(outfpath, 'wt') as outfile:

    # Write header
    outfile.write('ass_id\tdiff_large\tall_truncated\tpivotal_gene_seqID\tpivotal_gene_len\tmin_len\tmax_len\n')

    # Iterate over Assembly IDs
    for i, ass_id in enumerate(ass_ids):

        print(f'\rDoing {i+1}/{len(ass_ids)}: {ass_id}', end=' '*10)

        # Save some values to variables for speeeeeed
        curr_grpd_df = grpd_df[grpd_df['ass_id'] == ass_id]
        min_len = curr_grpd_df['min_len'].values[0]
        max_len = curr_grpd_df['max_len'].values[0]
        lendiff = curr_grpd_df['lendiff'].values[0]

        if lendiff > lendiff_threshold:
            # It genes differ in length greatly (> lendiff_threshold),
            #   we will search for pivotal genes (launch cmscan)

            # Select rows corresponding to corrent assembly
            curr_df = stats_df[stats_df['ass_id'] == ass_id]

            # ACCESSION.VERSION's
            accs = tuple(curr_df['acc'])

            # Select gene sequences of current genome
            select_seqs(accs, fasta_seqs_fpath, query_fasta_fpath)

            # Configure path to .tblout file
            tblout_fpath = os.path.join(tblout_dpath, f'{ass_id}.tblout')

            # Run cmscan, obtain `tblout_df`, which contains scores
            out_text, tblout_df = run_cmscan(cmscan_fpath, query_fasta_fpath, rfam_fpath, tblout_fpath, tblout_header)

            # We need seqIDs of extracted genes sequences
            seqIDs = parse_seqIDs(query_fasta_fpath)

            # Amend scores: penalize "local end alignments"
            amend_scores(seqIDs, out_text, tblout_df)

            # Update .tblout file
            tblout_df.to_csv(
                tblout_fpath,
                sep='\t',
                index=False,
                header=True,
                encoding='utf-8',
                na_rep='NA'
            )

            # Check if all genes are truncated
            all_truncated = 'no' not in set(tblout_df['trunc'])

            if all_truncated:
                # It all genes are truncated, we won't searhc for pivotal genes.
                # Write line indicating that genes are all aberrant
                #   and there is no need to search pivotal gene(s).
                outfile.write(f'{ass_id}\t1\t1\tNA\tNA\t{min_len}\t{max_len}\n')
            else:
                # It there is at least one non-truncated gene, we will select seqeunce
                #   with the hoghest score -- it will be pivotal.
                # There can be multiple pivotal genes, actually
                best_score = tblout_df['score'].max() - 1e-6
                best_gene_seqIDs = tuple(tblout_df[tblout_df['score'] >= best_score]['query_name'])

                # We want to see lengths of pivotal genes to verify ourselves
                pivotal_gene_len_dict = extract_pivotal_gene_lengths(query_fasta_fpath, best_gene_seqIDs)

                # Write result lines: one line per pivotal gene
                for seqID, length in pivotal_gene_len_dict.items():
                    outfile.write(f'{ass_id}\t1\t0\t{seqID}\t{length}\t{min_len}\t{max_len}\n')
                # end for
            # end if
        else:
            # Write line indicating that genes don't differ in lengths much
            #   and there is no need to search pivotal gene(s).
            outfile.write(f'{ass_id}\t0\tNA\tNA\tNA\t{min_len}\t{max_len}\n')
        # end if
    # end for
# end with

print('\nCompleted!')
print(outfpath)
