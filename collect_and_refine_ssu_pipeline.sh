
set -e

source corner_config.conf

for some_dir in "${WORKDIR}" "${LOGS_DIR}" "${GENOMES_GBK_DIR}"; do
  if [[ ! -d "${some_dir}" ]]; then
    mkdir -p "${some_dir}"
  fi
done


# == Translate Assembly UIDs to RefSeq GI numbers ==

ASS_ID_TO_GI_FPATH="${WORKDIR}/assembly_2_refseq.tsv"

# ./assembly2refseq_id.py \
#   --assm-id-file "${ASSEMBLY_IDS_FPATH}" \
#   --outfile "${ASS_ID_TO_GI_FPATH}"


# == Translate RefSeq GI numbers to corresponding ACCESSION.VERSION's and titles ==

GI_ACC_TITLES_FPATH="${WORKDIR}/archaea_refseq_accs.tsv"

# ./gis_to_accs.py \
#   --gi-file "${ASS_ID_TO_GI_FPATH}" \
#   --outfile "${GI_ACC_TITLES_FPATH}"


# == Merge Assembly IDs to ACCESSION.VERSION's and titles ==
# Moreover, this will remove "whole genome shotgun" sequences

ASS_ACC_MERGED_FPATH="${WORKDIR}/archaea_refseq_accs_merged.tsv"

# ./merge_assID2acc_and_remove_WGS.py \
#   --assm-2-gi-file "${ASS_ID_TO_GI_FPATH}" \
#   --gi-2-acc-file "${GI_ACC_TITLES_FPATH}" \
#   --outfile "${ASS_ACC_MERGED_FPATH}"


# == Download genomes ==

./download_genomes.py \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --outdir "${GENOMES_GBK_DIR}" \
  --log-file "${LOGS_DIR}/archaea_genome_download_log.log"


# == Extract 16S genes from downloaded genomes ==

ALL_GENES_FASTA="${WORKDIR}/gene_seqs/all_collected.fasta"
ALL_GENES_STATS="${WORKDIR}/gene_seqs/all_collected_stats.tsv"

./collect_16S/collect_16S.py \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --gbk-dir "${GENOMES_GBK_DIR}" \
  --out-fasta "${ALL_GENES_FASTA}" \
  --out-stats "${ALL_GENES_STATS}" \
  --cmsearch "${CMSEARCH_FOR_COLLECT_16S}" \
  --rfam-family-cm "${RFAM_FOR_COLLECT_16S}" \
  --seqkit "${SEQKIT}"
