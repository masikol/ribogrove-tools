#!/usr/bin/env bash

set -eu

# Get argumenets
CONF_FILE="$1"
ASM_ACC_TO_RM="$2"


# Check configuration file existance
if [[ ! -f "${CONF_FILE}" ]]; then
  echo -e "Error: configuration file ${CONF_FILE} does not exist!" >&2
  print_help
  exit 1
fi

# Load configuration
source "${CONF_FILE}"

if [[ "${DOMAIN}" != 'bacteria' && "${DOMAIN}" != 'archaea' ]]; then
  echo "Error: unrecognized domain: '${DOMAIN}'." >&2
  echo "It must be either 'bacteria' or 'archaea'." >&2
  exit 1
fi

if [[ -z "${ASM_ACC_TO_RM}" ]]; then
    echo "Error: please provide assembly accession to remove as the second argument" >&2
    exit 1
fi


GREP='grep'
if [[ ! -z $(which rg) ]]; then
    GREP='rg'
fi


TMP_FILE="${WORKDIR}/post_hoc_tmp.txt"

self_path=$(realpath "${0}")
self_dir_path=$(dirname ${self_path})
rm_from_json_gz_script="${self_dir_path}/remove_assembly_from_json.py"


echo "$(date) -- Starting"
echo "Working directory: '${WORKDIR}'"
echo "Assembly accession to remove: '${ASM_ACC_TO_RM}'"
echo '---------'


GENOMES_DATA_DIR="${WORKDIR}/genomes_data"
ASS_SUM_FILT_1="${GENOMES_DATA_DIR}/assembly_summary_filt1.txt.gz"
ASS_SUM_FINAL="${GENOMES_DATA_DIR}/assembly_summary_final.txt.gz"
REPLICON_MAP="${GENOMES_DATA_DIR}/replicon_map.tsv.gz"
gzipped_files_to_edit=(
    "${ASS_SUM_FILT_1}"
    "${ASS_SUM_FINAL}"
    "${REPLICON_MAP}"
)

CATEGORIES_DIR="${WORKDIR}/categories"
CATEGORIES_FILE="${CATEGORIES_DIR}/categories.tsv"
TAXONOMY_DIR="${WORKDIR}/taxonomy"
TAXONOMY_FILE="${TAXONOMY_DIR}/taxonomy.tsv"
GCNS_DIR="${WORKDIR}/GCNs"
PRIMERS_DIRPATH="${WORKDIR}/primers_coverage"
GENES_STATS_DIR="${WORKDIR}/gene_stats"
ABERRATIONS_AND_HETEROGENEITY_DIR="${WORKDIR}/aberrations_and_heterogeneity"
ABERRANT_SEQIDS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/aberrant_seqIDs.txt"
ENTROPY_FILE="${ABERRATIONS_AND_HETEROGENEITY_DIR}/entropy.tsv"
REPEATS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/repeats.tsv"
REPEAT_FAIL_SEQIDS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/repeats_fail_seqIDs.txt"
RIBOTYPER_FAIL_SEQIDS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/ribotyper_fail_seqIDs.txt"
RIBOTYPER_OUTDIR="${ABERRATIONS_AND_HETEROGENEITY_DIR}/ribotyper_out"


plain_files_to_edit=(
    "${CATEGORIES_FILE}"
    "${TAXONOMY_FILE}"
    ${GCNS_DIR}/*.tsv
    ${PRIMERS_DIRPATH}/*.tsv
    ${GENES_STATS_DIR}/*.tsv
    "${ABERRANT_SEQIDS_FPATH}"
    "${ENTROPY_FILE}"
    "${REPEATS_FPATH}"
    "${REPEAT_FAIL_SEQIDS_FPATH}"
    "${RIBOTYPER_FAIL_SEQIDS_FPATH}"
    "${ABERRATIONS_AND_HETEROGENEITY_DIR}/deletions.tsv"
    "${ABERRATIONS_AND_HETEROGENEITY_DIR}/insertions.tsv"
    "${ABERRATIONS_AND_HETEROGENEITY_DIR}/pident_pivotal_genes.tsv"
    "${ABERRATIONS_AND_HETEROGENEITY_DIR}/pivotal_genes.tsv"
    ${RIBOTYPER_OUTDIR}/*
)

GENES_DIR="${WORKDIR}/gene_seqs"
plain_fasta_files_to_edit=(
    ${GENES_DIR}/*.fasta
)

json_gz_files=(
    "${ABERRATIONS_AND_HETEROGENEITY_DIR}/per_base_entropy.json.gz"
)

echo "$(date) -- REMOVING FROM GZIPPED FILES"
for file in "${gzipped_files_to_edit[@]}"; do
    echo "$(date) -- removing from '${file}'"
    echo "  Lines before: $(zcat ${file} | wc -l)"
    zcat "${file}" \
        | "${GREP}" -Fv "${ASM_ACC_TO_RM}" \
        > "${TMP_FILE}"
    cat "${TMP_FILE}" \
        | gzip > "${file}"
    rm "${TMP_FILE}"
    echo "  Lines after: $(zcat ${file} | wc -l)"
done

echo -e "\n$(date) -- REMOVING FROM PLAIN TEXT FILES"
for file in "${plain_files_to_edit[@]}"; do
    echo "$(date) -- removing from '${file}'"
    echo "  Lines before: $(cat ${file} | wc -l)"
    cat "${file}" \
        | "${GREP}" -Fv "${ASM_ACC_TO_RM}" \
        > "${TMP_FILE}"
    cat "${TMP_FILE}" \
        > "${file}"
    rm "${TMP_FILE}"
    echo "  Lines after: $(cat ${file} | wc -l)"
done

echo -e "\n$(date) -- REMOVING FROM PLAIN TEXT FASTA FILES"
for file in "${plain_fasta_files_to_edit[@]}"; do
    echo "$(date) -- removing from '${file}'"
    echo '  Seq stats before:'
    seqkit stats "${file}"
    seqkit grep -nrvp "${ASM_ACC_TO_RM}" "${file}" \
        > "${TMP_FILE}"
    cat "${TMP_FILE}" > "${file}"
    rm "${TMP_FILE}"
    echo '  Seq stats after:'
    seqkit stats "${file}"
done

echo -e "\n$(date) -- REMOVING FROM JSON.GZ FILES"
for file in "${json_gz_files[@]}"; do
    echo "$(date) -- removing from '${file}'"
    python3 "${rm_from_json_gz_script}" \
        "${file}" \
        "${ASM_ACC_TO_RM}"
done


echo -e "\n--------------------"
echo 'Completed!'
echo 'Have fun and please come again!'
echo "Don't forget to add this assembly accession to the blacklist!"
exit 0
