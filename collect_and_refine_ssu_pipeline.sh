
set -e

source corner_config.conf


# == Translate Assembly UIDs to RefSeq GI numbers ==

ASS_ID_TO_GI_FPATH="${WORKDIR}/assembly_2_refseq.tsv"

# ./assembly2refseq_id.py \
#     --assm-id-fpath "${ASSEMBLY_IDS_FPATH}" \
#     --outfpath "${ASS_ID_TO_GI_FPATH}"


# == Translate RefSeq GI numbers to corresponding ACCESSION.VERSION's and titles ==

GI_ACC_TITLES_FPATH="${WORKDIR}/archaea_refseq_accs.tsv"

./gis_to_accs.py \
    --gi-fpath "${ASS_ID_TO_GI_FPATH}" \
    --outfpath "${GI_ACC_TITLES_FPATH}"
