
import os
import re
import sys

ASM_ACC_PATTERN = re.compile(
    r'([A-Z]{3})_([0-9]{3})([0-9]{3})([0-9]{3})\.([1-9][0-9]*)'
)


def get_asm_report_fpath(asm_accession, genomes_dirpath):
    return os.path.join(
        get_asm_data_dir_path(asm_accession, genomes_dirpath),
        '{}_assembly_report.txt'.format(asm_accession)
    )
# end def


def get_asm_data_dir_path(asm_accession, genomes_dirpath):
    # GCF_900012635.1
    # GCF             <- group(1)
    #     900         <- group(2)
    #        012      <- group(3)
    #           635   <- group(4)
    #               1 <- group(5)
    re_obj = re.match(ASM_ACC_PATTERN, asm_accession)
    if re_obj is None:
        sys.stderr.write(
            'Error: unexpected format of an assembly accession number: `{}`\n' \
                .format(asm_accession)
        )
        sys.exit(1)
    # end if
    dir_paths = [
        re_obj.group(1),
        re_obj.group(2),
        re_obj.group(3),
        re_obj.group(4),
        re_obj.group(5),
    ]
    return os.path.join(
        genomes_dirpath,
        os.path.join(*dir_paths)
    )
# end def


def get_genome_seqannot_fpath(asm_accession, genomes_dirpath):
    return os.path.join(
        get_asm_data_dir_path(asm_accession, genomes_dirpath),
        '{}.gbff.gz'.format(asm_accession)
    )
# end def


def primer_pair_key_2_outfpath(outdir_path, primer_pair_key):
    return os.path.join(outdir_path, f'{primer_pair_key}.tsv')
# end def
