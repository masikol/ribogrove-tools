
set -e

# Get work directory
WORKDIR="$1"
if [[ ! -d "${WORKDIR}" ]]; then
  echo -e "\nError: direcotry ${WORKDIR} does not exist!"
  exit 1
fi

# Get output directory
OUTDIR="$2"
if [[ ! -d "${OUTDIR}" ]]; then
  mkdir -v "${OUTDIR}"
fi

# Get release number
RELEASE_NUMBER="$3"
if [[ ! "${RELEASE_NUMBER}" =~ ^[0-9]+\.[0-9]+$ ]]; then
    echo "Release number is not of the required format. Allowed example: '2.208'"
    echo "This is the number specified by you: ${RELEASE_NUMBER}"
    exit 1
fi


# Dependencies

# https://github.com/shenwei356/csvtk
if [[ -z `which csvtk` ]]; then
    echo 'Cannot find csvtk program'
    echo 'Please, download it (https://github.com/shenwei356/csvtk) and add it to the PATH environment variable'
    exit 1
fi

# https://github.com/shenwei356/seqkit
if [[ -z `which seqkit` ]]; then
    echo 'Cannot find seqkit program'
    echo 'Please, download it (https://github.com/shenwei356/seqkit) and add it to the PATH environment variable'
    exit 1
fi

# Functions

function check_file() {
    fpath=$1
    if [[ ! -f ${fpath} ]]; then
        echo "Error: file '${fpath}' does not exist!"
        exit 1
    fi
}


# Variables for work

BACTERIA_DIR="${WORKDIR}/bacteria"
ARCHAEA_DIR="${WORKDIR}/archaea"

# Fasta

# Final fasta
BACTERIA_ANNOTATED_FASTA="${BACTERIA_DIR}/gene_seqs/bacteria_final_gene_seqs_annotated.fasta"
ARCHAEA_ANNOTATED_FASTA="${ARCHAEA_DIR}/gene_seqs/archaea_final_gene_seqs_annotated.fasta"
for f in "${BACTERIA_ANNOTATED_FASTA}" "${ARCHAEA_ANNOTATED_FASTA}"; do
    check_file "${f}"
done
RIBOGROVE_FASTA="${OUTDIR}/ribogrove_${RELEASE_NUMBER}_sequences.fasta.gz"


# == Metadata ==

METADATA_DIR="${OUTDIR}/metadata"
if [[ ! -d "${METADATA_DIR}" ]]; then
    mkdir -v "${METADATA_DIR}"
fi

# Discarded fasta
BACTERIA_DISCARDED_ANNOTATED_FASTA="${BACTERIA_DIR}/gene_seqs/bacteria_discarded_gene_seqs_annotated.fasta"
ARCHAEA_DISCARDED_ANNOTATED_FASTA="${ARCHAEA_DIR}/gene_seqs/archaea_discarded_gene_seqs_annotated.fasta"
for f in "${BACTERIA_DISCARDED_ANNOTATED_FASTA}" "${ARCHAEA_DISCARDED_ANNOTATED_FASTA}"; do
    check_file "${f}"
done
DISCARDED_FASTA="${METADATA_DIR}/discarded_${RELEASE_NUMBER}_sequences.fasta.gz"


# Index file
SCRIPT_DIR=`dirname $0`
INDEX_TEMPLATE="${SCRIPT_DIR}/INDEX_template.txt"
check_file "${INDEX_TEMPLATE}"
INDEX_FILE="${METADATA_DIR}/INDEX.txt"

# Source RefSeq genomes
BACTERIA_SOURCE_GENOMES="${BACTERIA_DIR}/bacteria_refseq_accs_merged.tsv"
ARCHAEA_SOURCE_GENOMES="${ARCHAEA_DIR}/archaea_refseq_accs_merged.tsv"
for f in "${BACTERIA_SOURCE_GENOMES}" "${ARCHAEA_SOURCE_GENOMES}"; do
    check_file "${f}"
done
RIBOGROVE_SOURCE_GENOMES="${METADATA_DIR}/source_RefSeq_genomes.tsv"

# Gene sequences statistics (final)
BACTERIA_PER_GENE_STATS="${BACTERIA_DIR}/bacteria_per_gene_stats.tsv"
ARCHAEA_PER_GENE_STATS="${ARCHAEA_DIR}/archaea_per_gene_stats.tsv"
for f in "${BACTERIA_PER_GENE_STATS}" "${ARCHAEA_PER_GENE_STATS}"; do
    check_file "${f}"
done
RIBOGROVE_PER_GENE_STATS="${METADATA_DIR}/gene_seqs_statistics.tsv"

# Gene sequences statistics (raw)
BACTERIA_DISCARDED_PER_GENE_STATS="${BACTERIA_DIR}/bacteria_discarded_per_gene_stats.tsv"
ARCHAEA_DISCARDED_PER_GENE_STATS="${ARCHAEA_DIR}/archaea_discarded_per_gene_stats.tsv"
for f in "${BACTERIA_DISCARDED_PER_GENE_STATS}" "${ARCHAEA_DISCARDED_PER_GENE_STATS}"; do
    check_file "${f}"
done
RIBOGROVE_DISCARDED_PER_GENE_STATS="${METADATA_DIR}/discarded_gene_seqs_statistics.tsv"

# Categories
BACTERIA_CATEGORIES="${BACTERIA_DIR}/categories/bacteria_categories.tsv"
ARCHAEA_CATEGORIES="${ARCHAEA_DIR}/categories/archaea_categories.tsv"
for f in "${BACTERIA_CATEGORIES}" "${ARCHAEA_CATEGORIES}"; do
    check_file "${f}"
done
RIBOGROVE_CATEGORIES="${METADATA_DIR}/categories.tsv"

# Taxonomy
BACTERIA_TAXONOMY="${BACTERIA_DIR}/taxonomy/bacteria_per_gene_taxonomy.tsv"
ARCHAEA_TAXONOMY="${ARCHAEA_DIR}/taxonomy/archaea_per_gene_taxonomy.tsv"
for f in "${BACTERIA_TAXONOMY}" "${ARCHAEA_TAXONOMY}"; do
    check_file "${f}"
done
RIBOGROVE_TAXONOMY="${METADATA_DIR}/taxonomy.tsv"

# Intragenic repeats
BACTERIA_REPEATS="${BACTERIA_DIR}/aberrations_and_heterogeneity/bacteria_repeats.tsv"
ARCHAEA_REPEATS="${ARCHAEA_DIR}/aberrations_and_heterogeneity/archaea_repeats.tsv"
for f in "${BACTERIA_REPEATS}" "${ARCHAEA_REPEATS}"; do
    check_file "${f}"
done
RIBOGROVE_REPEATS="${METADATA_DIR}/intragenic_repeats.tsv"

# Cmscan output table
BACTERIA_CMSCAN_TBLOUT="${BACTERIA_DIR}/aberrations_and_heterogeneity/cmscan_output_table.tblout"
ARCHAEA_CMSCAN_TBLOUT="${ARCHAEA_DIR}/aberrations_and_heterogeneity/cmscan_output_table.tblout"
for f in "${BACTERIA_CMSCAN_TBLOUT}" "${ARCHAEA_CMSCAN_TBLOUT}"; do
    check_file "${f}"
done
RIBOGROVE_TBLOUT="${METADATA_DIR}/cmscan_output_table.tblout"

# Entropy summary
BACTERIA_ENTROPY_SUMMARY="${BACTERIA_DIR}/aberrations_and_heterogeneity/bacteria_entropy.tsv"
ARCHAEA_ENTROPY_SUMMARY="${ARCHAEA_DIR}/aberrations_and_heterogeneity/archaea_entropy.tsv"
for f in "${BACTERIA_ENTROPY_SUMMARY}" "${ARCHAEA_ENTROPY_SUMMARY}"; do
    check_file "${f}"
done
RIBOGROVE_ENTROPY_SUMMARY="${METADATA_DIR}/entropy_summary.tsv"

# QIIME2-compatible taxonomy
QIIME2_TAXONOMY="${METADATA_DIR}/QIIME2-compatible-taxonomy.txt"
QIIME2_CONVERSION_SCRIPT='./make_qiime_taxonomy_file.py'
if [[ ! -f "${QIIME2_CONVERSION_SCRIPT}" ]]; then
    echo "Error: cannot find file '${QIIME2_CONVERSION_SCRIPT}'"
    exit 1
fi


# |======= PROCEED =======|

echo -e "\n|=== START PREPARING THE RELEASE ===|"

# ======= FASTA =======

# == Merge final annotated fasta files ==
# Sort them by sequences so that gzip will compress them even more

echo -n 'Final fasta...  '
if [[ -f "${RIBOGROVE_FASTA}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_FASTA}"
fi
cat "${BACTERIA_ANNOTATED_FASTA}" | seqkit sort -s | gzip > "${RIBOGROVE_FASTA}"
cat "${ARCHAEA_ANNOTATED_FASTA}" | seqkit sort -s | gzip >> "${RIBOGROVE_FASTA}"
echo "${RIBOGROVE_FASTA}"


# == Merge annotated fasta files of discarded sequences ==
# Sort them by sequences so that gzip will compress them even more

# Select discarded sequences
echo -n 'Discarded fasta...  '
if [[ -f "${DISCARDED_FASTA}" ]]; then
    # Empty the merged file
    echo -n '' > "${DISCARDED_FASTA}"
fi
cat "${BACTERIA_DISCARDED_ANNOTATED_FASTA}" \
    | seqkit sort -s | gzip > "${DISCARDED_FASTA}"
cat "${ARCHAEA_DISCARDED_ANNOTATED_FASTA}" \
    | seqkit sort -s | gzip >> "${DISCARDED_FASTA}"
echo "${DISCARDED_FASTA}"

# ======= METADATA =======

# Index file
echo -n 'Index file...  '
REPLACEMENT='RELEASE_NUMBER_REPLACE'
cat "${INDEX_TEMPLATE}" | sed "s|${REPLACEMENT}|${RELEASE_NUMBER}|g" > "${INDEX_FILE}"
echo "${INDEX_FILE}"

# Source RefSeq genomes
echo -n 'Source RefSeq genomes...  '
if [[ -f "${RIBOGROVE_SOURCE_GENOMES}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_SOURCE_GENOMES}"
fi
cat "${BACTERIA_SOURCE_GENOMES}" > "${RIBOGROVE_SOURCE_GENOMES}"
cat "${ARCHAEA_SOURCE_GENOMES}" | csvtk del-header -tT >> "${RIBOGROVE_SOURCE_GENOMES}"
echo "${RIBOGROVE_SOURCE_GENOMES}"

# Gene sequences statistics (final)
echo -n 'Gene sequences statistics (final)...  '
if [[ -f "${RIBOGROVE_PER_GENE_STATS}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_PER_GENE_STATS}"
fi
cat "${BACTERIA_PER_GENE_STATS}" > "${RIBOGROVE_PER_GENE_STATS}"
cat "${ARCHAEA_PER_GENE_STATS}" | csvtk del-header -tT >> "${RIBOGROVE_PER_GENE_STATS}"
echo "${RIBOGROVE_PER_GENE_STATS}"

# Gene sequences statistics (discarded)
echo -n 'Gene sequences statistics (discarded)...  '
if [[ -f "${RIBOGROVE_DISCARDED_PER_GENE_STATS}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_DISCARDED_PER_GENE_STATS}"
fi
cat "${BACTERIA_DISCARDED_PER_GENE_STATS}" > "${RIBOGROVE_DISCARDED_PER_GENE_STATS}"
cat "${ARCHAEA_DISCARDED_PER_GENE_STATS}" | csvtk del-header -tT >> "${RIBOGROVE_DISCARDED_PER_GENE_STATS}"
echo "${RIBOGROVE_DISCARDED_PER_GENE_STATS}"

# Categories
echo -n 'Categories...  '
if [[ -f "${RIBOGROVE_CATEGORIES}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_CATEGORIES}"
fi
cat "${BACTERIA_CATEGORIES}" > "${RIBOGROVE_CATEGORIES}"
cat "${ARCHAEA_CATEGORIES}" | csvtk del-header -tT >> "${RIBOGROVE_CATEGORIES}"
echo "${RIBOGROVE_CATEGORIES}"

# Taxonomy
echo -n 'Taxonomy...  '
if [[ -f "${RIBOGROVE_TAXONOMY}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_TAXONOMY}"
fi
cat "${BACTERIA_TAXONOMY}" > "${RIBOGROVE_TAXONOMY}"
cat "${ARCHAEA_TAXONOMY}" | csvtk del-header -tT >> "${RIBOGROVE_TAXONOMY}"
echo "${RIBOGROVE_TAXONOMY}"

# Intragenic repeats
echo -n 'Intragenic repeats...  '
if [[ -f "${RIBOGROVE_REPEATS}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_REPEATS}"
fi
cat "${BACTERIA_REPEATS}" > "${RIBOGROVE_REPEATS}"
cat "${ARCHAEA_REPEATS}" | csvtk del-header -tT >> "${RIBOGROVE_REPEATS}"
echo "${RIBOGROVE_REPEATS}"

# Cmscan output table
echo -n 'Cmscan output table...  '
if [[ -f "${RIBOGROVE_TBLOUT}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_TBLOUT}"
fi
cat "${BACTERIA_CMSCAN_TBLOUT}" > "${RIBOGROVE_TBLOUT}"
cat "${ARCHAEA_CMSCAN_TBLOUT}" | csvtk del-header -tT >> "${RIBOGROVE_TBLOUT}"
echo "${RIBOGROVE_TBLOUT}"

# Entropy summary
echo -n 'Entropy summary...  '
if [[ -f "${RIBOGROVE_ENTROPY_SUMMARY}" ]]; then
    # Empty the merged file
    echo -n '' > "${RIBOGROVE_ENTROPY_SUMMARY}"
fi
cat "${BACTERIA_ENTROPY_SUMMARY}" > "${RIBOGROVE_ENTROPY_SUMMARY}"
cat "${ARCHAEA_ENTROPY_SUMMARY}" | csvtk del-header -tT >> "${RIBOGROVE_ENTROPY_SUMMARY}"
echo "${RIBOGROVE_ENTROPY_SUMMARY}"

# QIIME2-compatible taxonomy
echo -n 'QIIME2-compatible taxonomy...  '
python3 "${QIIME2_CONVERSION_SCRIPT}" \
    --in-ribogrove-taxonomy "${RIBOGROVE_TAXONOMY}" \
    --out-qiime2-taxonomy "${QIIME2_TAXONOMY}"
echo "${QIIME2_TAXONOMY}"

# Zip the metadata
echo "Zipping the metadata directory: '${METADATA_DIR}/'"
WORKDIR=`pwd`
cd "${OUTDIR}"
METADATA_ARCHIVE="metadata_ribogrove_${RELEASE_NUMBER}.zip"
zip -r "${METADATA_ARCHIVE}" `basename "${METADATA_DIR}/"`
cd "${WORKDIR}"


echo -e "\n|=== END ===|"
echo "${OUTDIR}"
