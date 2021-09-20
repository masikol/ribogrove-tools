
set -e

source archaea_config.conf
# source bacteria_config.conf

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


# == Configure all paths to intermediate and result files ==

ASS_ID_TO_GI_FPATH="${WORKDIR}/${PREFIX}_assembly_2_refseq.tsv"

GI_ACC_TITLES_FPATH="${WORKDIR}/${PREFIX}_refseq_accs.tsv"

ASS_ACC_MERGED_FPATH="${WORKDIR}/${PREFIX}_refseq_accs_merged.tsv"

ALL_GENES_FASTA="${GENES_DIR}/${PREFIX}_all_collected.fasta"
ALL_GENES_STATS="${GENES_DIR}/${PREFIX}_all_collected_stats.tsv"

PER_GENOME_CAT_FPATH="${CATEGORIES_DIR}/${PREFIX}_per_genome_categories.tsv"
PER_GENE_CAT_FPATH="${CATEGORIES_DIR}/${PREFIX}_per_gene_categories.tsv"
SEQTECH_LOGFILE="${LOGS_DIR}/${PREFIX}_seqtech_log.log"

PER_GENOME_TAXID_FPATH="${TAXONOMY_DIR}/${PREFIX}_per_genome_taxIDs.tsv"

PER_GENE_TAXID_FPATH="${TAXONOMY_DIR}/${PREFIX}_per_gene_taxIDs.tsv"

PER_GENOME_TAXONOMY_FPATH="${TAXONOMY_DIR}/${PREFIX}_per_genome_taxonomy.tsv"
PER_GENE_TAXONOMY_FPATH="${TAXONOMY_DIR}/${PREFIX}_per_gene_taxonomy.tsv"

NO_NN_FASTA_FPATH="${GENES_DIR}/${PREFIX}_gene_seqs_no_NN.fasta"
NO_NN_STATS_FPATH="${GENES_DIR}/${PREFIX}_gene_stats_no_NN.tsv"
NN_FASTA_FPATH="${GENES_DIR}/${PREFIX}_NN_gene_seqs.fasta"

REPEATS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/${PREFIX}_repeats.tsv"

PIVOTAL_GENES_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/${PREFIX}_pivotal_genes.tsv"

ABERRANT_SEQIDS_FPATH="${ABERRATIONS_AND_HETEROGENEITY_DIR}/${PREFIX}_aberrant_seqIDs.txt"
NON_ABERRANT_GENES_FASTA="${GENES_DIR}/${PREFIX}_non_aberrant_gene_seqs.fasta"
NON_ABERRANT_GENES_STATS="${GENES_DIR}/${PREFIX}_non_aberrant_gene_stats.tsv"

PURE_GENES_FASTA="${GENES_DIR}/${PREFIX}_pure_gene_seqs.fasta"
PURE_GENES_STATS="${GENES_DIR}/${PREFIX}_pure_gene_stats.tsv"


# == Translate Assembly UIDs to RefSeq GI numbers ==

./assembly2refseq_id.py \
  --assm-id-file "${ASSEMBLY_IDS_FPATH}" \
  --outfile "${ASS_ID_TO_GI_FPATH}"


# == Translate RefSeq GI numbers to corresponding ACCESSION.VERSION's and titles ==

./gis_to_accs.py \
  --gi-file "${ASS_ID_TO_GI_FPATH}" \
  --outfile "${GI_ACC_TITLES_FPATH}"


# == Merge Assembly IDs to ACCESSION.VERSION's and titles ==
# Moreover, this will remove "whole genome shotgun" sequences

./merge_assID2acc_and_remove_WGS.py \
  --assm-2-gi-file "${ASS_ID_TO_GI_FPATH}" \
  --gi-2-acc-file "${GI_ACC_TITLES_FPATH}" \
  --outfile "${ASS_ACC_MERGED_FPATH}"


# == Download genomes ==

./download_genomes.py \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --outdir "${GENOMES_GBK_DIR}" \
  --log-file "${LOGS_DIR}/archaea_genome_download_log.log"


# == Extract 16S genes from downloaded genomes ==

./collect_16S/collect_16S.py \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --gbk-dir "${GENOMES_GBK_DIR}" \
  --out-fasta "${ALL_GENES_FASTA}" \
  --out-stats "${ALL_GENES_STATS}" \
  --cmsearch "${CMSEARCH_FOR_COLLECT_16S}" \
  --rfam-family-cm "${RFAM_FOR_COLLECT_16S}" \
  --seqkit "${SEQKIT}"


# == Assign categories to downloaded genomes ==

./assign_genome_categories/assign_genome_categories.py \
  --all-fasta-file "${ALL_GENES_FASTA}" \
  --all-stats-file "${ALL_GENES_STATS}" \
  --gbk-dir "${GENOMES_GBK_DIR}" \
  --per-genome-outfile "${PER_GENOME_CAT_FPATH}" \
  --per-gene-outfile "${PER_GENE_CAT_FPATH}" \
  --seqtech-logfile "${SEQTECH_LOGFILE}" \
  --seqkit "${SEQKIT}"


# == Get taxIDs for our genomes ==

./get_taxIDs.py \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --all-fasta-file "${ALL_GENES_FASTA}" \
  --per-genome-outfile "${PER_GENOME_TAXID_FPATH}" \
  --seqkit "${SEQKIT}"


# == Map seqIDs to taxIDs ==

./pergenome_2_pergene_taxIDs.py \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --all-fasta-file "${ALL_GENES_FASTA}" \
  --per-genome-taxID-file "${PER_GENOME_TAXID_FPATH}" \
  --per-gene-outfile "${PER_GENE_TAXID_FPATH}" \


# == Map our Aseembly IDs (and seqIDs) to full taxonomy using our taxIDs ==

./add_taxonomy_names.py \
  --per-genome-taxid-file "${PER_GENOME_TAXID_FPATH}" \
  --per-gene-taxid-file "${PER_GENE_TAXID_FPATH}" \
  --ranked-lineage "${RANKEDLINEAGE_FPATH}" \
  --per-genome-outfile "${PER_GENOME_TAXONOMY_FPATH}" \
  --per-gene-outfile "${PER_GENE_TAXONOMY_FPATH}" \
  --seqkit "${SEQKIT}"


# == Drop genes containing at leats 2 N's in a row ==

./drop_NN.py \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --all-fasta-file "${ALL_GENES_FASTA}" \
  --out-fasta-file "${NO_NN_FASTA_FPATH}" \
  --out-stats-file "${NO_NN_STATS_FPATH}" \
  --NN-outfile "${NN_FASTA_FPATH}"


# == Find repeats in genes sequences ==

./find_repeats.py \
  --no-NN-fasta-file "${NO_NN_FASTA_FPATH}" \
  --conserved-regions-fasta "${CONSERVED_REGIONS_FASTA}" \
  --outfile "${REPEATS_FPATH}"


# == Find pivotal genes ==

./find_pivotal_genes.py \
  --fasta-seqs-file "${NO_NN_FASTA_FPATH}" \
  --genes-stats-file "${NO_NN_STATS_FPATH}" \
  --outfile "${PIVOTAL_GENES_FPATH}" \
  --tblout-dir "${TBLOUT_DIR}" \
  --cmscan "${CMSCAN_FOR_PIVOTAL_GENES}" \
  --cmpress "${CMPRESS_FOR_PIVOTAL_GENES}" \
  --rfam-family-cm "${RFAM_FOR_PIVOTAL_GENES}" \
  --lendiff-threshold 5


# == Find aberrant genes and record long indels ==

./find_aberrant_genes.py \
  --fasta-seqs-file "${NO_NN_FASTA_FPATH}" \
  --genes-stats-file "${NO_NN_STATS_FPATH}" \
  --pivotal-genes-file "${PIVOTAL_GENES_FPATH}" \
  --conserved-regions-fasta "${CONSERVED_REGIONS_FASTA}" \
  --outdir "${ABERRATIONS_AND_HETEROGENEITY_DIR}" \
  --muscle "${MUSCLE}" \
  --indel-len-threshold 10


# == Drop aberarant genes ==

./drop_aberrant_genes.py \
  --input-fasta-file "${NO_NN_FASTA_FPATH}" \
  --aberrant-genes-file "${ABERRANT_SEQIDS_FPATH}" \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --categories-file "${PER_GENE_CAT_FPATH}" \
  --output-fasta "${NON_ABERRANT_GENES_FASTA}" \
  --output-stats "${NON_ABERRANT_GENES_STATS}" \
  --seqkit "${SEQKIT}"


# == Drop long repeats ==

./drop_repeats.py \
  --input-fasta-file "${NON_ABERRANT_GENES_FASTA}" \
  --genes-stats-file "${NON_ABERRANT_GENES_STATS}" \
  --assm-acc-file "${ASS_ACC_MERGED_FPATH}" \
  --repeats-file "${REPEATS_FPATH}" \
  --categories-file "${PER_GENE_CAT_FPATH}" \
  --out-fasta "${PURE_GENES_FASTA}" \
  --out-stats "${PURE_GENES_STATS}"