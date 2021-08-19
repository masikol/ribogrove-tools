
set -e

source corner_config.conf


# == Translate Assembly UIDs to RefSeq GI numbers ==

ASS_ID_TO_GI_FPATH="${WORKDIR}/assembly_2_refseq.tsv"

# ./assembly2refseq_id.py \
#     --assm-id-fpath "${ASSEMBLY_IDS_FPATH}" \
#     --outfpath "${ASS_ID_TO_GI_FPATH}"


# == Translate RefSeq GI numbers to corresponding ACCESSION.VERSION's and titles ==

GI_ACC_TITLES_FPATH="${WORKDIR}/archaea_refseq_accs.tsv"

# ./gis_to_accs.py \
#     --gi-fpath "${ASS_ID_TO_GI_FPATH}" \
#     --outfpath "${GI_ACC_TITLES_FPATH}"


# == Merge Assembly IDs to ACCESSION.VERSION's and titles ==
# Moreover, this will remove "whole genome shotgun" sequences

ASS_ACC_MERGED_FPATH="${WORKDIR}/archaea_refseq_accs_merged.tsv"

./merge_assID2acc_and_remove_WGS.py \
    --assm-2-gi-fpath "${ASS_ID_TO_GI_FPATH}" \
    --gi-2-acc-fpath "${GI_ACC_TITLES_FPATH}" \
    --outfpath "${ASS_ACC_MERGED_FPATH}"
