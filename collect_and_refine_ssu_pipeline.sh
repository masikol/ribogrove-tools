
source corner_config.conf

ASS_ID_TO_GI_FPATH="${WORKDIR}/assembly_2_refseq.tsv"

./assembly2refseq_id.py \
    --assm-id-fpath "${ASSEMBLY_IDS_FPATH}" \
    --outfpath "${ASS_ID_TO_GI_FPATH}"
