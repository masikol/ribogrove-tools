
import os


def get_asm_report_fpath(asm_accession, genomes_dirpath):
    return os.path.join(
        genomes_dirpath,
        asm_accession,
        '{}_assembly_report.txt'.format(asm_accession)
    )
# end def


def get_genome_seqannot_fpath(asm_accession, genomes_dirpath):
    return os.path.join(
        genomes_dirpath,
        asm_accession,
        '{}.gbff.gz'.format(asm_accession)
    )
# end def


def primer_pair_key_2_outfpath(outdir_path, primer_pair_key):
    return os.path.join(outdir_path, f'{primer_pair_key}.tsv')
# end def
