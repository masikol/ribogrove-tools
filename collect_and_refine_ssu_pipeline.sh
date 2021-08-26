
set -e

# source corner_config.conf
source bacteria_config.conf

# Directories for different sorts of data
LOGS_DIR="${WORKDIR}/logs"
CATEGORIES_DIR="${WORKDIR}/categories"
TAXONOMY_DIR="${WORKDIR}/taxonomy"
GENES_DIR="${WORKDIR}/gene_seqs"
ABERRATIONS_AND_HETEROGENEITY_DIR="${WORKDIR}/aberrations_and_heterogeneity"
TBLOUT_DIR="${WORKDIR}/cmscan_tblout"

for some_dir in "${WORKDIR}" "${LOGS_DIR}" "${CATEGORIES_DIR}" \
                "${TAXONOMY_DIR}" "${GENOMES_GBK_DIR}" "${GENES_DIR}" \
                "${ABERRATIONS_AND_HETEROGENEITY_DIR}" "${TBLOUT_DIR}";
do
  if [[ ! -d "${some_dir}" ]]; then
    mkdir -pv "${some_dir}"
  fi
done


# == Translate Assembly UIDs to RefSeq GI numbers ==

ASS_ID_TO_GI_FPATH="${WORKDIR}/assembly_2_refseq.tsv"

# ./assembly2refseq_id.py \
#   --assm-id-file "${ASSEMBLY_IDS_FPATH}" \
#   --outfile "${ASS_ID_TO_GI_FPATH}"


# == Translate RefSeq GI numbers to corresponding ACCESSION.VERSION's and titles ==

# GI_ACC_TITLES_FPATH="${WORKDIR}/archaea_refseq_accs.tsv"
GI_ACC_TITLES_FPATH="${WORKDIR}/bacteria_ass_refseq_accs.tsv"

# ./gis_to_accs.py \
#   --gi-file "${ASS_ID_TO_GI_FPATH}" \
#   --outfile "${GI_ACC_TITLES_FPATH}"


# == Merge Assembly IDs to ACCESSION.VERSION's and titles ==
# Moreover, this will remove "whole genome shotgun" sequences

# ASS_ACC_MERGED_FPATH="${WORKDIR}/archaea_refseq_accs_merged.tsv"
ASS_ACC_MERGED_FPATH="${WORKDIR}/bacteria_ass_refseq_accs_merged.tsv"

# ./merge_assID2acc_and_remove_WGS.py \
#   --assm-2-gi-file "${ASS_ID_TO_GI_FPATH}" \
#   --gi-2-acc-file "${GI_ACC_TITLES_FPATH}" \
#   --outfile "${ASS_ACC_MERGED_FPATH}"


# == Download genomes ==

# ./download_genomes.py \
#   --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
#   --outdir "${GENOMES_GBK_DIR}" \
#   --log-file "${LOGS_DIR}/archaea_genome_download_log.log"


# == Extract 16S genes from downloaded genomes ==

ALL_GENES_FASTA="${GENES_DIR}/all_collected.fasta"
ALL_GENES_STATS="${GENES_DIR}/all_collected_stats.tsv"

# ./collect_16S/collect_16S.py \
#   --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
#   --gbk-dir "${GENOMES_GBK_DIR}" \
#   --out-fasta "${ALL_GENES_FASTA}" \
#   --out-stats "${ALL_GENES_STATS}" \
#   --cmsearch "${CMSEARCH_FOR_COLLECT_16S}" \
#   --rfam-family-cm "${RFAM_FOR_COLLECT_16S}" \
#   --seqkit "${SEQKIT}"


# == Assign categories to downloaded genomes ==

# PER_GENOME_CAT_FPATH="${CATEGORIES_DIR}/archaea_per_genome_categories.tsv"
# PER_GENE_CAT_FPATH="${CATEGORIES_DIR}/archaea_per_gene_categories.tsv"
# SEQTECH_LOGFILE="${LOGS_DIR}/archaea_seqtech_log.log"
PER_GENOME_CAT_FPATH="${CATEGORIES_DIR}/bacteria_per_genome_categories.tsv"
PER_GENE_CAT_FPATH="${CATEGORIES_DIR}/bacteria_per_gene_categories.tsv"
SEQTECH_LOGFILE="${LOGS_DIR}/bacteria_seqtech_log.log"

# ./assign_genome_categories/assign_genome_categories.py \
#   --all-fasta-file "${ALL_GENES_FASTA}" \
#   --all-stats-file "${ALL_GENES_STATS}" \
#   --gbk-dir "${GENOMES_GBK_DIR}" \
#   --per-genome-outfile "${PER_GENOME_CAT_FPATH}" \
#   --per-gene-outfile "${PER_GENE_CAT_FPATH}" \
#   --seqtech-logfile "${SEQTECH_LOGFILE}" \
#   --seqkit "${SEQKIT}"


# == Get taxIDs for our genomes ==

# PER_GENOME_TAXID_FPATH="${TAXONOMY_DIR}/archaea_per_genome_taxIDs.tsv"
# PER_GENE_TAXID_FPATH="${TAXONOMY_DIR}/archaea_per_gene_taxIDs.tsv"
PER_GENOME_TAXID_FPATH="${TAXONOMY_DIR}/bacteria_per_genome_taxIDs.tsv"
PER_GENE_TAXID_FPATH="${TAXONOMY_DIR}/bacteria_per_gene_taxIDs.tsv"

# ./get_taxIDs.py \
#   --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
#   --all-fasta-file "${ALL_GENES_FASTA}" \
#   --per-genome-outfile "${PER_GENOME_TAXID_FPATH}" \
#   --per-gene-outfile "${PER_GENE_TAXID_FPATH}" \
#   --seqkit "${SEQKIT}"


# == Map our Aseembly IDs (and seqIDs) to full taxonomy using our taxIDs ==

# PER_GENOME_TAXONOMY_FPATH="${TAXONOMY_DIR}/archaea_per_genome_taxonomy.tsv"
# PER_GENE_TAXONOMY_FPATH="${TAXONOMY_DIR}/archaea_per_gene_taxonomy.tsv"
PER_GENOME_TAXONOMY_FPATH="${TAXONOMY_DIR}/bacteria_per_genome_taxonomy.tsv"
PER_GENE_TAXONOMY_FPATH="${TAXONOMY_DIR}/bacteria_per_gene_taxonomy.tsv"

./add_taxonomy_names.py \
  --per-genome-taxid-file "${PER_GENOME_TAXID_FPATH}" \
  --per-gene-taxid-file "${PER_GENE_TAXID_FPATH}" \
  --ranked-lineage "${RANKEDLINEAGE_FPATH}" \
  --per-genome-outfile "${PER_GENOME_TAXONOMY_FPATH}" \
  --per-gene-outfile "${PER_GENE_TAXONOMY_FPATH}" \
  --seqkit "${SEQKIT}"


# == Drop genes containing at leats 2 N's in a row ==

NO_NN_FASTA_FPATH="${GENES_DIR}/gene_seqs_no_NN.fasta"
NO_NN_STATS_FPATH="${GENES_DIR}/gene_stats_no_NN.tsv"
NN_FASTA_FPATH="${GENES_DIR}/NN_gene_seqs.fasta"

# ./drop_NN.py \
#   --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
#   --all-fasta-file "${ALL_GENES_FASTA}" \
#   --out-fasta-file "${NO_NN_FASTA_FPATH}" \
#   --out-stats-file "${NO_NN_STATS_FPATH}" \
#   --NN-outfile "${NN_FASTA_FPATH}"


# == Find repeats in genes sequences ==

REPEATS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/repeats.tsv"

# ./find_repeats.py \
#   --no-NN-fasta-file "${NO_NN_FASTA_FPATH}" \
#   --conserved-regions-fasta "${CONSERVED_REGIONS_FASTA}" \
#   --outfile "${REPEATS_FPATH}"


# == Find pivotal genes ==

PIVOTAL_GENES_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/pivotal_genes.tsv"

# ./find_pivotal_genes.py \
#   --fasta-seqs-file "${NO_NN_FASTA_FPATH}" \
#   --genes-stats-file "${NO_NN_STATS_FPATH}" \
#   --outfile "${PIVOTAL_GENES_FPATH}" \
#   --tblout-dir "${TBLOUT_DIR}" \
#   --cmscan "${CMSCAN_FOR_PIVOTAL_GENES}" \
#   --cmpress "${CMPRESS_FOR_PIVOTAL_GENES}" \
#   --rfam-family-cm "${RFAM_FOR_PIVOTAL_GENES}" \
#   --lendiff-threshold 5



# == Find aberrant genes and record long indels ==

# ./find_aberrant_genes.py \
#   --fasta-seqs-file "${NO_NN_FASTA_FPATH}" \
#   --genes-stats-file "${NO_NN_STATS_FPATH}" \
#   --pivotal-genes-file "${PIVOTAL_GENES_FPATH}" \
#   --conserved-regions-fasta "${CONSERVED_REGIONS_FASTA}" \
#   --outdir "${ABERRATIONS_AND_HETEROGENEITY_DIR}" \
#   --muscle "${MUSCLE}" \
#   --indel-len-threshold 10


# == Drop aberarant genes ==

ABERRANT_SEQIDS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/aberrant_seqIDs.txt"
NON_ABERRANT_GENES_FASTA="${GENES_DIR}/non_aberrant_gene_seqs.fasta"
NON_ABERRANT_GENES_STATS="${GENES_DIR}/non_aberrant_genes_stats.tsv"

# ./drop_aberrant_genes.py \
#   --input-fasta-file "${NO_NN_FASTA_FPATH}" \
#   --aberrant-genes-file "${ABERRANT_SEQIDS_FPATH}" \
#   --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
#   --categories-file "${PER_GENE_CAT_FPATH}" \
#   --output-fasta "${NON_ABERRANT_GENES_FASTA}" \
#   --output-stats "${NON_ABERRANT_GENES_STATS}" \
#   --seqkit "${SEQKIT}"
