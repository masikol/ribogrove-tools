set -e


function print_help {
  echo 'Usage:' >&2
  echo "  $0 <WORKDIR_ROOT> <RELEASE_OUTDIR> <RELEASE_NUMBER>" >&2
  echo 'Example:' >&2
  echo "  $0 RiboGrove_9.215_workdir RiboGrove_9.215_release 9.215" >&2
}

if [[ ${#} -eq 0 || $1 == '-h' || $1 == '-help' || $1 == '--help' ]]; then
  print_help
  exit 1
fi


# Get work directory
WORKDIR="$1"
if [[ ! -d "${WORKDIR}" ]]; then
  echo -e "\nError: directory ${WORKDIR} does not exist!" >&2
  print_help
  exit 1
fi

# Get output directory
OUTDIR="$2"
if [[ ! -d "${OUTDIR}" ]]; then
  mkdir -vp "${OUTDIR}"
  if [[ $? != 0 ]]; then
    echo "Error: cannot reate output dir '${OUTDIR}'" >&2
    exit 1
  fi
fi

# Get release number
RELEASE_NUMBER="$3"
if [[ ! "${RELEASE_NUMBER}" =~ ^[0-9]+\.[0-9]+$ ]]; then
  echo "Release number is not of the required format. Allowed example: '2.208'" >&2
  echo "And here is the release you've specified: '${RELEASE_NUMBER}"'' >&2
  print_help
  exit 1
fi


# Dependencies

# https://github.com/shenwei356/csvtk
if [[ -z `which csvtk` ]]; then
  echo 'Cannot find csvtk program' >&2
  echo 'Please, download it (https://github.com/shenwei356/csvtk)' >&2
  echo '  and add it to the PATH environment variable' >&2
  exit 1
fi

# https://github.com/shenwei356/seqkit
if [[ -z `which seqkit` ]]; then
  echo 'Cannot find seqkit program' >&2
  echo 'Please, download it (https://github.com/shenwei356/seqkit)' >&2
  echo '  and add it to the PATH environment variable' >&2
  exit 1
fi


# Functions

function check_file() {
  fpath=$1
  if [[ ! -f ${fpath} ]]; then
    echo "Error: file '${fpath}' does not exist!" >&2
    exit 1
  fi
}


# |== Set variables -- path to input and output files ==|

BACTERIA_DIR="${WORKDIR}/bacteria"
ARCHAEA_DIR="${WORKDIR}/archaea"

# Fasta

# Final fasta
BACTERIA_ANNOTATED_FASTA="${BACTERIA_DIR}/gene_seqs/final_gene_seqs_annotated.fasta"
ARCHAEA_ANNOTATED_FASTA="${ARCHAEA_DIR}/gene_seqs/final_gene_seqs_annotated.fasta"
for f in "${BACTERIA_ANNOTATED_FASTA}" "${ARCHAEA_ANNOTATED_FASTA}"; do
  check_file "${f}"
done
RIBOGROVE_FASTA="${OUTDIR}/ribogrove_${RELEASE_NUMBER}_sequences.fasta.gz"


# == Metadata ==

METADATA_DIR="${OUTDIR}/metadata"
if [[ ! -d "${METADATA_DIR}" ]]; then
  mkdir -vp "${METADATA_DIR}"
  if [[ $? != 0 ]]; then
  echo "Error: cannot reate metadata dir '${METADATA_DIR}'" >&2
  exit 1
  fi
fi

# Discarded fasta
BACTERIA_DISCARDED_ANNOTATED_FASTA="${BACTERIA_DIR}/gene_seqs/discarded_gene_seqs_annotated.fasta"
ARCHAEA_DISCARDED_ANNOTATED_FASTA="${ARCHAEA_DIR}/gene_seqs/discarded_gene_seqs_annotated.fasta"
for f in "${BACTERIA_DISCARDED_ANNOTATED_FASTA}" "${ARCHAEA_DISCARDED_ANNOTATED_FASTA}"; do
  check_file "${f}"
done
DISCARDED_FASTA="${METADATA_DIR}/discarded_${RELEASE_NUMBER}_sequences.fasta.gz"


# Index file
SCRIPT_DIR=`dirname $0`
INDEX_TEMPLATE="${SCRIPT_DIR}/templates/INDEX_template.txt"
check_file "${INDEX_TEMPLATE}"
INDEX_FILE="${METADATA_DIR}/INDEX.txt"

# Source RefSeq genomes
BACTERIA_ASM_SUM="${BACTERIA_DIR}/genomes_data/assembly_summary_final.txt.gz"
ARCHAEA_ASM_SUM="${ARCHAEA_DIR}/genomes_data/assembly_summary_final.txt.gz"
for f in "${BACTERIA_ASM_SUM}" "${ARCHAEA_ASM_SUM}"; do
  check_file "${f}"
done
RIBOGROVE_SOURCE_GENOMES="${METADATA_DIR}/source_RefSeq_genomes.tsv"

# Gene sequences statistics (final)
BACTERIA_PER_GENE_STATS="${BACTERIA_DIR}/gene_stats/per_gene_stats.tsv"
ARCHAEA_PER_GENE_STATS="${ARCHAEA_DIR}/gene_stats/per_gene_stats.tsv"
for f in "${BACTERIA_PER_GENE_STATS}" "${ARCHAEA_PER_GENE_STATS}"; do
    check_file "${f}"
done
RIBOGROVE_PER_GENE_STATS="${METADATA_DIR}/gene_seqs_statistics.tsv"

# Gene sequences statistics (raw)
BACTERIA_DISCARDED_PER_GENE_STATS="${BACTERIA_DIR}/gene_stats/discarded_per_gene_stats.tsv"
ARCHAEA_DISCARDED_PER_GENE_STATS="${ARCHAEA_DIR}/gene_stats/discarded_per_gene_stats.tsv"
for f in "${BACTERIA_DISCARDED_PER_GENE_STATS}" "${ARCHAEA_DISCARDED_PER_GENE_STATS}"; do
    check_file "${f}"
done
RIBOGROVE_DISCARDED_PER_GENE_STATS="${METADATA_DIR}/discarded_gene_seqs_statistics.tsv"

# Categories
BACTERIA_CATEGORIES="${BACTERIA_DIR}/categories/categories.tsv"
ARCHAEA_CATEGORIES="${ARCHAEA_DIR}/categories/categories.tsv"
for f in "${BACTERIA_CATEGORIES}" "${ARCHAEA_CATEGORIES}"; do
  check_file "${f}"
done
RIBOGROVE_CATEGORIES="${METADATA_DIR}/categories.tsv"

# Taxonomy
BACTERIA_TAXONOMY="${BACTERIA_DIR}/taxonomy/taxonomy.tsv"
ARCHAEA_TAXONOMY="${ARCHAEA_DIR}/taxonomy/taxonomy.tsv"
for f in "${BACTERIA_TAXONOMY}" "${ARCHAEA_TAXONOMY}"; do
  check_file "${f}"
done
RIBOGROVE_TAXONOMY="${METADATA_DIR}/taxonomy.tsv"

# Intragenic repeats
BACTERIA_REPEATS="${BACTERIA_DIR}/aberrations_and_heterogeneity/repeats.tsv"
ARCHAEA_REPEATS="${ARCHAEA_DIR}/aberrations_and_heterogeneity/repeats.tsv"
for f in "${BACTERIA_REPEATS}" "${ARCHAEA_REPEATS}"; do
  check_file "${f}"
done
RIBOGROVE_REPEATS="${METADATA_DIR}/intragenic_repeats.tsv"
REPEAT_LEN_THRESHOLD=25

# Entropy summary
BACTERIA_ENTROPY_SUMMARY="${BACTERIA_DIR}/aberrations_and_heterogeneity/entropy.tsv"
ARCHAEA_ENTROPY_SUMMARY="${ARCHAEA_DIR}/aberrations_and_heterogeneity/entropy.tsv"
for f in "${BACTERIA_ENTROPY_SUMMARY}" "${ARCHAEA_ENTROPY_SUMMARY}"; do
  check_file "${f}"
done
RIBOGROVE_ENTROPY_SUMMARY="${METADATA_DIR}/entropy_summary.tsv"


# |======= PROCEED =======|

echo -e "\n|=== START PREPARING THE RELEASE ===|"

# ------- FINAL SEQUENCES -------

# == Merge final annotated fasta files ==
# Sort them by sequences so that gzip will compress them even more

echo -n 'Final fasta...  '
if [[ -f "${RIBOGROVE_FASTA}" ]]; then
  # Empty the merged file
  echo -n '' > "${RIBOGROVE_FASTA}"
fi
cat "${BACTERIA_ANNOTATED_FASTA}" "${ARCHAEA_ANNOTATED_FASTA}" \
  | seqkit sort -s \
  | gzip \
  > "${RIBOGROVE_FASTA}"
echo "${RIBOGROVE_FASTA}"


# ------- METADATA -------

# Index file
echo -n 'Index file...  '
REPLACEMENT='RELEASE_NUMBER_REPLACE'
cat "${INDEX_TEMPLATE}" | sed "s|${REPLACEMENT}|${RELEASE_NUMBER}|g" > "${INDEX_FILE}"
echo "${INDEX_FILE}"

# == Merge annotated fasta files of discarded sequences ==
# Sort them by sequences so that gzip will compress them even more
echo -n 'Discarded fasta...  '
if [[ -f "${DISCARDED_FASTA}" ]]; then
  # Empty the merged file
  echo -n '' > "${DISCARDED_FASTA}"
fi
cat "${BACTERIA_DISCARDED_ANNOTATED_FASTA}" "${ARCHAEA_DISCARDED_ANNOTATED_FASTA}" \
  | seqkit sort -s \
  | gzip \
  > "${DISCARDED_FASTA}"
echo "${DISCARDED_FASTA}"

# Source RefSeq genomes
echo -n 'Source RefSeq genomes...  '
if [[ -f "${RIBOGROVE_SOURCE_GENOMES}" ]]; then
  # Empty the merged file
  echo -n '' > "${RIBOGROVE_SOURCE_GENOMES}"
fi
columns='asm_acc,taxid,species_taxid,organism_name,infraspecific_name,assembly_level,genome_rep,seq_rel_date,asm_name,submitter'
zcat "${BACTERIA_ASM_SUM}" \
  | csvtk cut -tTf "${columns}" \
  > "${RIBOGROVE_SOURCE_GENOMES}"
cat "${ARCHAEA_ASM_SUM}" \
  | csvtk cut -tTf "${columns}" \
  | csvtk del-header -tT \
  >> "${RIBOGROVE_SOURCE_GENOMES}"
echo "${RIBOGROVE_SOURCE_GENOMES}"

# Gene sequences statistics (final)
echo -n 'Gene sequences statistics (final)...  '
if [[ -f "${RIBOGROVE_PER_GENE_STATS}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_PER_GENE_STATS}"
fi
cat "${BACTERIA_PER_GENE_STATS}" > "${RIBOGROVE_PER_GENE_STATS}"
cat "${ARCHAEA_PER_GENE_STATS}" \
  | csvtk del-header -tT \
  >> "${RIBOGROVE_PER_GENE_STATS}"
echo "${RIBOGROVE_PER_GENE_STATS}"

# Gene sequences statistics (discarded)
echo -n 'Gene sequences statistics (discarded)...  '
if [[ -f "${RIBOGROVE_DISCARDED_PER_GENE_STATS}" ]]; then
  # Empty the merged file
  echo -n '' > "${RIBOGROVE_DISCARDED_PER_GENE_STATS}"
fi
cat "${BACTERIA_DISCARDED_PER_GENE_STATS}" > "${RIBOGROVE_DISCARDED_PER_GENE_STATS}"
cat "${ARCHAEA_DISCARDED_PER_GENE_STATS}" \
  | csvtk del-header -tT \
  >> "${RIBOGROVE_DISCARDED_PER_GENE_STATS}"
echo "${RIBOGROVE_DISCARDED_PER_GENE_STATS}"

# Categories
echo -n 'Categories...  '
if [[ -f "${RIBOGROVE_CATEGORIES}" ]]; then
  # Empty the merged file
  echo -n '' > "${RIBOGROVE_CATEGORIES}"
fi
cat "${BACTERIA_CATEGORIES}" > "${RIBOGROVE_CATEGORIES}"
cat "${ARCHAEA_CATEGORIES}" \
  | csvtk del-header -tT \
  >> "${RIBOGROVE_CATEGORIES}"
echo "${RIBOGROVE_CATEGORIES}"

# Taxonomy
echo -n 'Taxonomy...  '
if [[ -f "${RIBOGROVE_TAXONOMY}" ]]; then
  # Empty the merged file
  echo -n '' > "${RIBOGROVE_TAXONOMY}"
fi
cat "${BACTERIA_TAXONOMY}" > "${RIBOGROVE_TAXONOMY}"
cat "${ARCHAEA_TAXONOMY}" \
  | csvtk del-header -tT \
  >> "${RIBOGROVE_TAXONOMY}"
echo "${RIBOGROVE_TAXONOMY}"

# Intragenic repeats
echo -n 'Intragenic repeats...  '
if [[ -f "${RIBOGROVE_REPEATS}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_REPEATS}"
fi
FILTER_STR="rep_len>${REPEAT_LEN_THRESHOLD}"
cat "${BACTERIA_REPEATS}" \
  | csvtk filter -tTf "${FILTER_STR}" \
  > "${RIBOGROVE_REPEATS}"
cat "${ARCHAEA_REPEATS}" \
  | csvtk filter -tTf "${FILTER_STR}" \
  | csvtk del-header -tT \
  >> "${RIBOGROVE_REPEATS}"
echo "${RIBOGROVE_REPEATS}"

# Entropy summary
echo -n 'Entropy summary...  '
if [[ -f "${RIBOGROVE_ENTROPY_SUMMARY}" ]]; then
  # Empty the merged file
  echo -n '' > "${RIBOGROVE_ENTROPY_SUMMARY}"
fi
cat "${BACTERIA_ENTROPY_SUMMARY}" > "${RIBOGROVE_ENTROPY_SUMMARY}"
cat "${ARCHAEA_ENTROPY_SUMMARY}" \
  | csvtk del-header -tT \
  >> "${RIBOGROVE_ENTROPY_SUMMARY}"
echo "${RIBOGROVE_ENTROPY_SUMMARY}"

# Zip the metadata
echo "Zipping the metadata directory: '${METADATA_DIR}/'"
cd "${OUTDIR}"
METADATA_ARCHIVE="metadata_ribogrove_${RELEASE_NUMBER}.zip"
zip -r "${METADATA_ARCHIVE}" `basename "${METADATA_DIR}"`


echo -e "\n|=== END ===|"
echo "${OUTDIR}"
