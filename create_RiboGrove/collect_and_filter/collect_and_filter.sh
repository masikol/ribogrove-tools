set -e

function print_help {
  echo 'Usage:' >&2
  echo "  bash ${0} <CONFIG_FILE> [-ta]" >&2
  echo '-t: Test mode. The script will not use entire RefSeq' >&2
  echo '  but just some genomes from the file' >&2
  echo '  scripts/data/test/test_bacteria_assembly_summary.txt.gz or' >&2
  echo '  scripts/data/test/test_archaea_assembly_summary.txt.gz' >&2
  echo '-a: RefSeq catalog file is already filtered.' >&2
  echo '  If you already have file RefSeq-release216_filtered.catalog.gz' >&2
  echo '  created from file RefSeq-release216.catalog.gz,' >&2
  echo '  the script will not filter it again and thus save some time.' >&2
  echo 'Example:' >&2
  echo "  bash ${0} config/archaea_example.conf -a" >&2
}


if [[ ${#} -eq 0 || $1 == '-h' || $1 == '-help' || $1 == '--help' ]]; then
   print_help
   exit 1
fi


# Check configuration file existance
CONF_FILE="$1"
if [[ ! -f "${CONF_FILE}" ]]; then
  echo -e "Error: configuration file ${CONF_FILE} does not exist!" >&2
  print_help
  exit 1
fi


# Parse options
TEST_MODE=false
REFSEQ_CATALOG_ALREADY_FILTERED=false

shift
while getopts ":ta" opt; do
  case "${opt}" in
  t)
      TEST_MODE=true
      echo 'INFO: Running in test mode' >&2
      ;;
  a)
      REFSEQ_CATALOG_ALREADY_FILTERED=true
      echo 'INFO: Assuming that RefSeq catalog is already filtered.' >&2
      ;;
  esac
done

# Load configuration
source "${CONF_FILE}"

if [[ "${DOMAIN}" != 'bacteria' && "${DOMAIN}" != 'archaea' ]]; then
  echo "Error: unrecognized domain: '${DOMAIN}'." >&2
  echo "It must be either 'bacteria' or 'archaea'." >&2
  exit 1
fi


# Get directory where scripts and their data are located
script_abspath=`realpath $0`
pipeline_sh_dir=`dirname "${script_abspath}"`
SCRIPTS_DIR="${pipeline_sh_dir}/scripts"
SCRIPTS_DATA_DIR="${SCRIPTS_DIR}/data"


# Directories for different sorts of data
GENOMES_DATA_DIR="${WORKDIR}/genomes_data"
GENES_DIR="${WORKDIR}/gene_seqs"
GENES_STATS_DIR="${WORKDIR}/gene_stats"
CATEGORIES_DIR="${WORKDIR}/categories"
TAXONOMY_DIR="${WORKDIR}/taxonomy"
LOGS_DIR="${WORKDIR}/logs"
ABERRATIONS_AND_HETEROGENEITY_DIR="${WORKDIR}/aberrations_and_heterogeneity"
GCNS_DIR="${WORKDIR}/GCNs"

dirs_to_create=(
  "${WORKDIR}" "${LOGS_DIR}" "${CATEGORIES_DIR}" "${GCNS_DIR}" \
  "${TAXONOMY_DIR}" "${GENOMES_DATA_DIR}" "${GENOMES_GBK_DIR}" \
  "${GENES_DIR}" "${GENES_STATS_DIR}" "${ABERRATIONS_AND_HETEROGENEITY_DIR}"
)
if [[ "${CALC_PRIMERS_COVERAGE}" == true ]]; then
  PRIMERS_DIRPATH="${WORKDIR}/primers_coverage"
  dirs_to_create+=( "${PRIMERS_DIRPATH}" )
fi
for some_dir in ${dirs_to_create[@]}; do
  mkdir -pv "${some_dir}"
done



# |=== Configure all paths to intermediate and result files ===|

FILTERED_REFSEQ_CATALOG_FILE="${REFSEQ_CATALOG_FILE/.catalog.gz/_filtered.catalog.gz}"

if [[ "${TEST_MODE}" == false ]]; then
  ASS_SUM_LINK="https://ftp.ncbi.nlm.nih.gov/genomes/refseq/${DOMAIN}/assembly_summary.txt"
  ASS_SUM="${GENOMES_DATA_DIR}/assembly_summary.gz"
else
  ASS_SUM="${SCRIPTS_DATA_DIR}/test/test_${DOMAIN}_assembly_summary.txt.gz"
fi

ASS_SUM_FILT_1="${GENOMES_DATA_DIR}/assembly_summary_filt1.txt.gz"
ASS_SUM_FINAL="${GENOMES_DATA_DIR}/assembly_summary_final.txt.gz"
REPLICON_MAP="${GENOMES_DATA_DIR}/replicon_map.tsv.gz"

ASM_BLACKLIST_FPATH="${SCRIPTS_DATA_DIR}/ad_hoc/assembly_blacklist.tsv"
BLACKLIST_SEQIDS_FILE="${SCRIPTS_DATA_DIR}/ad_hoc/blacklist_seqIDs.tsv"
WHITELIST_SEQIDS_FILE="${SCRIPTS_DATA_DIR}/ad_hoc/whitelist_seqIDs.tsv"

ALL_GENES_FASTA="${GENES_DIR}/all_collected.fasta"
ALL_GENES_STATS="${GENES_STATS_DIR}/all_collected_stats.tsv"

TAXONOMY_FILE="${TAXONOMY_DIR}/taxonomy.tsv"

CATEGORIES_FILE="${CATEGORIES_DIR}/categories.tsv"

RIBOTYPER_OUTDIR="${ABERRATIONS_AND_HETEROGENEITY_DIR}/ribotyper_out"
RIBOTYPER_SHORT_OUT_TSV="${RIBOTYPER_OUTDIR}/ribotyper_out.ribotyper.short.out.tsv"
RIBOTYPER_LONG_OUT_TSV="${RIBOTYPER_OUTDIR}/ribotyper_out.ribotyper.long.out.tsv"
RIBOTYPER_FAIL_SEQIDS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/ribotyper_fail_seqIDs.txt"
RIBOTYPER_ACCEPT_FILE="${SCRIPTS_DATA_DIR}/${DOMAIN}.accept"

REPEATS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/repeats.tsv"
REPEAT_FAIL_SEQIDS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/repeats_fail_seqIDs.txt"

CMSCAN_TBLOUT_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/cmscan_output_table.tblout"

ABERRANT_SEQIDS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/aberrant_seqIDs.txt"

FINAL_GENES_FASTA="${GENES_DIR}/final_gene_seqs.fasta"

ANNOTATED_RESULT_FASTA="${GENES_DIR}/final_gene_seqs_annotated.fasta"
FINAL_SEQIDS="${GENES_DIR}/final_seqIDs.txt"
DISCARDED_FASTA="${GENES_DIR}/discarded_gene_seqs.fasta"
ANNOTATED_DISCARDED_FASTA="${GENES_DIR}/discarded_gene_seqs_annotated.fasta"

COUNT_BASES_TABLE="${GENES_STATS_DIR}/base_counts.tsv"
DISCARDED_COUNT_BASES_TABLE="${GENES_STATS_DIR}/discarded_base_counts.tsv"

ENTROPY_FILE="${ABERRATIONS_AND_HETEROGENEITY_DIR}/entropy.tsv"

if [[ ! -z "${PREV_WORKDIR}" ]]; then
  CACHE_MODE=true
  prev_aberr_dir="${PREV_WORKDIR}/aberrations_and_heterogeneity"
  PREV_ASM_SUM_FINAL="${PREV_WORKDIR}/genomes_data/assembly_summary_final.txt.gz"
  PREV_REPLICON_MAP="${PREV_WORKDIR}/genomes_data/replicon_map.tsv.gz"
  PREV_ASM_ACCS_NNN="${PREV_WORKDIR}/genomes_data/asm_accs_NNN.txt.gz"
  PREV_RIBOTYPER_SHORT_OUT_TSV="${prev_aberr_dir}/ribotyper_out/ribotyper_out.ribotyper.short.out.tsv"
  PREV_RIBOTYPER_LONG_OUT_TSV="${prev_aberr_dir}/ribotyper_out/ribotyper_out.ribotyper.long.out.tsv"
  PREV_TBLOUT_FILE="${prev_aberr_dir}/cmscan_output_table.tblout"
  PREV_ALL_GENES_FASTA="${PREV_WORKDIR}/gene_seqs/all_collected.fasta"
  PREV_ALL_GENES_STATS="${PREV_WORKDIR}/gene_stats/all_collected_stats.tsv"
  PREV_FINAL_GENES_FASTA="${PREV_WORKDIR}/gene_seqs/final_gene_seqs_annotated.fasta"
  PREV_PERBASE_ENTROPY_FILE="${prev_aberr_dir}/per_base_entropy.tsv.gz"
  PREV_PRIMERS_DIR="${PREV_WORKDIR}/primers_coverage"
  PREV_ABERRANT_SEQIDS="${prev_aberr_dir}/aberrant_seqIDs.txt"
else
  CACHE_MODE=false
fi


# |=== Proceed ===|

# == Filter RefSeq .catalog file ==

if [[ "${REFSEQ_CATALOG_ALREADY_FILTERED}" == false ]]; then
  python3 "${SCRIPTS_DIR}/filter_refseq_catalog.py" \
    --raw-refseq-catalog "${REFSEQ_CATALOG_FILE}" \
    --outfile "${FILTERED_REFSEQ_CATALOG_FILE}"
fi


# == Download assembly_summary.txt ==

if [[ "${TEST_MODE}" == false ]]; then
  curl \
    "${ASS_SUM_LINK}" \
    | gzip > "${ASS_SUM}"
fi


# == Filter assembly summary ==

python3 "${SCRIPTS_DIR}/filter_asm_summary_step1.py" \
  --in-ass-sum "${ASS_SUM}" \
  --refseq-catalog "${FILTERED_REFSEQ_CATALOG_FILE}" \
  --acc-blacklist "${ASM_BLACKLIST_FPATH}" \
  --out-ass-sum "${ASS_SUM_FILT_1}"


# == Download genomes ==

python3 "${SCRIPTS_DIR}/download_genomes.py" \
  --asm-sum "${ASS_SUM_FILT_1}" \
  --outdir "${GENOMES_GBK_DIR}" \
  --log-file "${LOGS_DIR}/download_genomes.log"


# == Make replicon map ==

if [[ "${CACHE_MODE}" == true ]]; then
  python3 "${SCRIPTS_DIR}/make_replicon_map.py" \
    --asm-sum "${ASS_SUM_FILT_1}" \
    --genomes-dir "${GENOMES_GBK_DIR}" \
    --prev-replicon-map "${PREV_REPLICON_MAP}" \
    --out "${REPLICON_MAP}"
else
  python3 "${SCRIPTS_DIR}/make_replicon_map.py" \
    --asm-sum "${ASS_SUM_FILT_1}" \
    --genomes-dir "${GENOMES_GBK_DIR}" \
    --out "${REPLICON_MAP}"
fi


# == Make final Assembly summary file ==
if [[ "${CACHE_MODE}" == true ]]; then
  python3 "${SCRIPTS_DIR}/filter_asm_summary_step2.py" \
    --in-asm-sum "${ASS_SUM_FILT_1}" \
    --replicon-map "${REPLICON_MAP}" \
    --refseq-catalog "${FILTERED_REFSEQ_CATALOG_FILE}" \
    --genomes-dir "${GENOMES_GBK_DIR}" \
    --prev-asm-sum-final "${PREV_ASM_SUM_FINAL}" \
    --prev-NNN-asm-accs "${PREV_ASM_ACCS_NNN}" \
    --out-asm-sum "${ASS_SUM_FINAL}"
else
  python3 "${SCRIPTS_DIR}/filter_asm_summary_step2.py" \
    --in-asm-sum "${ASS_SUM_FILT_1}" \
    --replicon-map "${REPLICON_MAP}" \
    --refseq-catalog "${FILTERED_REFSEQ_CATALOG_FILE}" \
    --genomes-dir "${GENOMES_GBK_DIR}" \
    --out-asm-sum "${ASS_SUM_FINAL}"
fi


# == Make taxonomy ==

python3 "${SCRIPTS_DIR}/make_taxonomy.py" \
    --asm-sum "${ASS_SUM_FINAL}" \
    --ranked-lineage "${RANKEDLINEAGE_FPATH}" \
    --out "${TAXONOMY_FILE}"


# == Extract 16S genes from downloaded genomes ==
if [[ "${CACHE_MODE}" == true ]]; then
  python3 "${SCRIPTS_DIR}/extract_16S.py" \
    --asm-sum "${ASS_SUM_FINAL}" \
    --genomes-dir "${GENOMES_GBK_DIR}" \
    --out-fasta "${ALL_GENES_FASTA}" \
    --out-stats "${ALL_GENES_STATS}" \
    --cmsearch "${CMSEARCH_FOR_EXTRACT_16S}" \
    --rfam-family-cm "${RFAM_FAMILY_FOR_EXTRACT_16S}" \
    --prev-all-genes-fasta "${PREV_ALL_GENES_FASTA}" \
    --prev-all-genes-stats "${PREV_ALL_GENES_STATS}" \
    --seqkit "${SEQKIT}"
else
  python3 "${SCRIPTS_DIR}/extract_16S.py" \
    --asm-sum "${ASS_SUM_FINAL}" \
    --genomes-dir "${GENOMES_GBK_DIR}" \
    --out-fasta "${ALL_GENES_FASTA}" \
    --out-stats "${ALL_GENES_STATS}" \
    --cmsearch "${CMSEARCH_FOR_EXTRACT_16S}" \
    --rfam-family-cm "${RFAM_FAMILY_FOR_EXTRACT_16S}" \
    --seqkit "${SEQKIT}"
fi


# == Assign categories to downloaded genomes ==

python3 "${SCRIPTS_DIR}/assign_genome_categories.py" \
  --all-fasta-file "${ALL_GENES_FASTA}" \
  --all-stats-file "${ALL_GENES_STATS}" \
  --genomes-dir "${GENOMES_GBK_DIR}" \
  --outfile "${CATEGORIES_FILE}" \
  --seqkit "${SEQKIT}"


# == Find inappropriate sequences with ribotyper ==

if [[ "${CACHE_MODE}" == true ]]; then
  python3 "${SCRIPTS_DIR}/check_seqs_with_ribotyper.py" \
    --in-fasta-file "${ALL_GENES_FASTA}" \
    --outdir "${RIBOTYPER_OUTDIR}" \
    --ribotyper "${RIBOTYPER}" \
    --ribotyper-threads "${RIBOTYPER_THREADS}" \
    --prev-short-out-tsv "${PREV_RIBOTYPER_SHORT_OUT_TSV}" \
    --prev-long-out-tsv "${PREV_RIBOTYPER_LONG_OUT_TSV}" \
    --acccept-file "${RIBOTYPER_ACCEPT_FILE}"
else
  python3 "${SCRIPTS_DIR}/check_seqs_with_ribotyper.py" \
    --in-fasta-file "${ALL_GENES_FASTA}" \
    --outdir "${RIBOTYPER_OUTDIR}" \
    --ribotyper "${RIBOTYPER}" \
    --ribotyper-threads "${RIBOTYPER_THREADS}" \
    --acccept-file "${RIBOTYPER_ACCEPT_FILE}"
fi

python3 "${SCRIPTS_DIR}/find_ribotyper_fail_seqs.py" \
  --in-short-out-tsv "${RIBOTYPER_SHORT_OUT_TSV}" \
  --out-fail-file "${RIBOTYPER_FAIL_SEQIDS_FPATH}"


# == Find aberrant genes and record long indels ==

if [[ "${CACHE_MODE}" == true ]]; then
  python3 "${SCRIPTS_DIR}/find_aberrant_genes.py" \
    --fasta-seqs-file "${ALL_GENES_FASTA}" \
    --ribotyper-fail-seqIDs "${RIBOTYPER_FAIL_SEQIDS_FPATH}" \
    --in-asm-sum "${ASS_SUM_FINAL}" \
    --ribotyper-long-out-tsv "${RIBOTYPER_LONG_OUT_TSV}" \
    --prev-final-fasta "${PREV_FINAL_GENES_FASTA}" \
    --prev-aberrant-seqIDs "${PREV_ABERRANT_SEQIDS}" \
    --outdir "${ABERRATIONS_AND_HETEROGENEITY_DIR}" \
    --mafft "${MAFFT}" \
    --deletion-len-threshold "${DELETION_LEN_THRESHOLD}"
else
  python3 "${SCRIPTS_DIR}/find_aberrant_genes.py" \
    --fasta-seqs-file "${ALL_GENES_FASTA}" \
    --ribotyper-fail-seqIDs "${RIBOTYPER_FAIL_SEQIDS_FPATH}" \
    --in-asm-sum "${ASS_SUM_FINAL}" \
    --ribotyper-long-out-tsv "${RIBOTYPER_LONG_OUT_TSV}" \
    --outdir "${ABERRATIONS_AND_HETEROGENEITY_DIR}" \
    --mafft "${MAFFT}" \
    --deletion-len-threshold "${DELETION_LEN_THRESHOLD}"
fi

# == Find repeats in genes sequences ==

python3 "${SCRIPTS_DIR}/find_repeats.py" \
  --in-fasta-file "${ALL_GENES_FASTA}" \
  --ribotyper-fail-seqIDs "${RIBOTYPER_FAIL_SEQIDS_FPATH}" \
  --aberrant-seqIDs "${ABERRANT_SEQIDS_FPATH}" \
  --repeat-len-threshold "${REPEAT_LEN_THRESHOLD}" \
  --out-fail-file "${REPEAT_FAIL_SEQIDS_FPATH}" \
  --out-repeats-log "${REPEATS_FPATH}"



# == Drop sequences which didn't pass filters ==
python3 "${SCRIPTS_DIR}/make_final_seqs.py" \
  --all-seqs-file "${ALL_GENES_FASTA}" \
  --ribotyper-fail-seqIDs "${RIBOTYPER_FAIL_SEQIDS_FPATH}" \
  --aberrant-seqIDs "${ABERRANT_SEQIDS_FPATH}" \
  --repeats-fail-seqIDs "${REPEAT_FAIL_SEQIDS_FPATH}" \
  --blacklist-seqIDs "${BLACKLIST_SEQIDS_FILE}" \
  --whitelist-seqIDs "${WHITELIST_SEQIDS_FILE}" \
  --replicon-map "${REPLICON_MAP}" \
  --out-fasta-file "${FINAL_GENES_FASTA}"



# == Annotate sequences: add taxonomy and categories to their headers ==
python3 "${SCRIPTS_DIR}/annotate_seq_names.py" \
  --fasta-seqs-file "${FINAL_GENES_FASTA}" \
  --taxonomy-file "${TAXONOMY_FILE}" \
  --categories-file "${CATEGORIES_FILE}" \
  --outfile "${ANNOTATED_RESULT_FASTA}"


# == Make result sequences pretty: 70 bp per line ==

tmp_fasta="${WORKDIR}/tmp.fasta"
cat "${ANNOTATED_RESULT_FASTA}" | "${SEQKIT}" seq -uw 70 > "${tmp_fasta}"
cat "${tmp_fasta}" > "${ANNOTATED_RESULT_FASTA}"
rm -v "${tmp_fasta}"


# == Make a fasta file of discarded sequences ==

cat "${ANNOTATED_RESULT_FASTA}" \
  | seqkit seq -ni > "${FINAL_SEQIDS}"
cat "${ALL_GENES_FASTA}" \
  | seqkit grep -vf "${FINAL_SEQIDS}" \
  | seqkit sort -s > "${DISCARDED_FASTA}"

if [[ -f "${FINAL_SEQIDS}" ]]; then
  rm -v "${FINAL_SEQIDS}"
fi


# == Annotate discarded sequences: add taxonomy and categories to their headers ==

python3 "${SCRIPTS_DIR}/annotate_seq_names.py" \
  --fasta-seqs-file "${DISCARDED_FASTA}" \
  --taxonomy-file "${TAXONOMY_FILE}" \
  --categories-file "${CATEGORIES_FILE}" \
  --outfile "${ANNOTATED_DISCARDED_FASTA}"


# == Make result sequences pretty: 70 bp per line ==

tmp_fasta="${WORKDIR}/tmp.fasta"
cat "${ANNOTATED_DISCARDED_FASTA}" | "${SEQKIT}" seq -uw 70 > "${tmp_fasta}"
cat "${tmp_fasta}" > "${ANNOTATED_DISCARDED_FASTA}"
rm -v "${tmp_fasta}"


# == Count bases (final) ==

python3 "${SCRIPTS_DIR}/count_bases.py" \
  --input-fasta "${ANNOTATED_RESULT_FASTA}" \
  --outfile "${COUNT_BASES_TABLE}"

# == Count bases for discarded sequences ==

python3 "${SCRIPTS_DIR}/count_bases.py" \
  --input-fasta "${ANNOTATED_DISCARDED_FASTA}" \
  --outfile "${DISCARDED_COUNT_BASES_TABLE}"


# == Calculate entropy -- intragenomic variability ==

if [[ "${CACHE_MODE}" == true ]]; then
  python3 "${SCRIPTS_DIR}/calculate_entropy.py" \
    --fasta-seqs-file "${ANNOTATED_RESULT_FASTA}" \
    --categories-file "${CATEGORIES_FILE}" \
    --outfile "${ENTROPY_FILE}" \
    --prev-per-base-entropy-file "${PREV_PERBASE_ENTROPY_FILE}" \
    --mafft "${MAFFT}"
else
  python3 "${SCRIPTS_DIR}/calculate_entropy.py" \
    --fasta-seqs-file "${ANNOTATED_RESULT_FASTA}" \
    --categories-file "${CATEGORIES_FILE}" \
    --outfile "${ENTROPY_FILE}" \
    --mafft "${MAFFT}"
fi


# == Calculate PCR primer coverage ==

python3 "${SCRIPTS_DIR}/check_primers_mfeprimer.py" \
  --fasta-seqs-file "${ALL_GENES_FASTA}" \
  --outdir "${PRIMERS_DIRPATH}" \
  --mfeprimer "${MFEPRIMER}" \
  --mfe-tmp-dir "${MFEPRIMER_TMP_DIR}" \
  --threads "${MFEPRIMER_THREADS}"



# == Calculate Gene Copy Numbers ==

if [[ "${CALC_PRIMERS_COVERAGE}" == true ]]; then
  python3 "${SCRIPTS_DIR}/calculate_GCNs.py" \
    --final-base-counts "${COUNT_BASES_TABLE}" \
    --primers-dir "${PRIMERS_DIRPATH}" \
    --outdir "${GCNS_DIR}"
else
  python3 "${SCRIPTS_DIR}/calculate_GCNs.py" \
    --final-base-counts "${DISCARDED_COUNT_BASES_TABLE}" \
    --outdir "${GCNS_DIR}"
fi
