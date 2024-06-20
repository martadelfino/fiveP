################################################################################
###### Script for obtaining gene paralogues ####################################
################################################################################

library('tidyverse')

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("biomaRt")
library('biomaRt')

#source('./src/functions.R')

# Install an old version of dplyr for biomaRt to work.
install.packages('https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.3.4.tar.gz', repos = NULL)


# Reading the protein coding genes file ----------------------------------------

#hgnc_ensembl <- read_delim('./data/raw/protein_coding_genes.txt', '\t') %>%
 # dplyr::select(hgnc_id, ensembl_gene_id)
#print(hgnc_ensembl)

rdata_hgnc_ensembl <- load('data/protein_coding_genes.RData')
hgnc_ensembl <- get("protein_coding_genes") %>%
  dplyr::select(hgnc_id, ensembl_gene_id)

ensembl_id_vector <- hgnc_ensembl %>% pull(ensembl_gene_id)

# Querying BioBart for gene paralogues -----------------------------------------
# threshold is 30%

#ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
#attributes = listAttributes(ensembl)
#print(attributes$page %>% unique())
#print(attributes %>% dplyr::filter(grepl('paralog', name)))

#name                                              description     page
#1                  hsapiens_paralog_ensembl_gene                           Human paralogue gene stable ID homologs
#2          hsapiens_paralog_associated_gene_name                     Human paralogue associated gene name homologs
#3               hsapiens_paralog_ensembl_peptide                 Human paralogue protein or transcript ID homologs
#4                    hsapiens_paralog_chromosome                 Human paralogue chromosome/scaffold name homologs
#5                   hsapiens_paralog_chrom_start           Human paralogue chromosome/scaffold start (bp) homologs
#6                     hsapiens_paralog_chrom_end             Human paralogue chromosome/scaffold end (bp) homologs
#7  hsapiens_paralog_canonical_transcript_protein                 Paralogue query protein or transcript ID homologs
#8                       hsapiens_paralog_subtype                Paralogue last common ancestor with Human homologs
#9                hsapiens_paralog_orthology_type                            Human paralogue homology type homologs
#10                      hsapiens_paralog_perc_id Paralogue %id. target Human gene identical to query gene homologs
#11                   hsapiens_paralog_perc_id_r1 Paralogue %id. query gene identical to target Human gene homologs


# I got these from a file, don't remember which one. but there are more?
#hsapiens_paralog_ensembl_gene                  Human Paralog Ensembl Gene ID
#hsapiens_paralog_canonical_transcript_protein  Canonical Protein or Transcript ID
#hsapiens_paralog_ensembl_peptide               Human Paralog Ensembl Protein ID
#hsapiens_paralog_chromosome                    Human Paralog Chromosome Name
#hsapiens_paralog_chrom_start                   Human Paralog Chr Start (bp)
#hsapiens_paralog_chrom_end                     Human Paralog Chr End (bp)
#hsapiens_paralog_orthology_type                Homology Type
#hsapiens_paralog_subtype                       Ancestor
#hsapiens_paralog_paralogy_confidence           Paralogy confidence [0 low, 1 high]
#hsapiens_paralog_perc_id                       % Identity with respect to query gene
#hsapiens_paralog_perc_id_r1                    % Identity with respect to Human gene
#hsapiens_paralog_dn                            dN
#hsapiens_paralog_ds                            dS


human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
paralogues <- biomaRt::getBM(attributes = c("ensembl_gene_id",
                                "hsapiens_paralog_ensembl_gene",
                                "hsapiens_paralog_orthology_type",
                                "hsapiens_paralog_perc_id",
                                "hsapiens_paralog_perc_id_r1",
                                "version"),
                 filters = "ensembl_gene_id",
                 values = ensembl_id_vector,
                 mart = human)
 head(paralogues)


# Saving the raw paralogues file -----------------------------------------------

save(paralogues, file = "data/biomart_paralogues.RData")

#save_data_with_datetime('./data/raw/', 'biomart_paralogues', paralogues)


# Cleaning the paralogues file -------------------------------------------------

# Map back the input gene ensembl IDs back to HGNC IDs
paralogues_cleaned <- paralogues %>%
  dplyr::left_join(hgnc_ensembl, by = 'ensembl_gene_id') %>%
  dplyr::select(hgnc_id, ensembl_gene_id, hsapiens_paralog_orthology_type,
                hsapiens_paralog_ensembl_gene, hsapiens_paralog_perc_id,
                hsapiens_paralog_perc_id_r1) %>%
  dplyr::rename(gene1_hgnc_id = hgnc_id) %>%
  dplyr::rename(gene1_ensembl_gene_id = ensembl_gene_id)
head(paralogues_cleaned)

# Now map the paralog ensembl IDs to HGNC IDs
paralogues_cleaned_with_paralog_hgnc <- paralogues_cleaned %>%
  left_join(hgnc_ensembl, join_by(hsapiens_paralog_ensembl_gene == ensembl_gene_id),
            relationship = "many-to-many") %>%
  dplyr::rename(paralog_hgnc_id = hgnc_id) %>%
  dplyr::rename(paralog_ensembl_gene_id = hsapiens_paralog_ensembl_gene) %>%
  dplyr::rename(paralog_perc_id = hsapiens_paralog_perc_id) %>%
  dplyr::rename(paralog_perc_id_r1 = hsapiens_paralog_perc_id_r1)
head(paralogues_cleaned_with_paralog_hgnc)

# Cleaning it up
paralogues_cleaned_reorg <- paralogues_cleaned_with_paralog_hgnc %>%
  dplyr::select(gene1_hgnc_id, gene1_ensembl_gene_id, paralog_hgnc_id,
                paralog_ensembl_gene_id, hsapiens_paralog_orthology_type,
                paralog_perc_id, paralog_perc_id_r1) %>%
  dplyr::arrange(gene1_hgnc_id)
head(paralogues_cleaned_reorg)


# Taking the mean of the two paralog percentages
paralogues_cleaned_reorg <- paralogues_cleaned_reorg %>%
  dplyr::mutate(mean_paralog_perc = (paralog_perc_id + paralog_perc_id_r1) / 2) %>%
  dplyr::mutate(max_paralog_perc = pmax(paralog_perc_id, paralog_perc_id_r1))
head(paralogues_cleaned_reorg)


# Saving the cleaned paralogues file -----------------------------------------------

#save_data_with_datetime('./data/processed/', 'biomart_paralogues_cleaned',
 #                       paralogues_cleaned_reorg)

save(paralogues_cleaned_reorg, file = "data/biomart_paralogues_cleaned.RData")





# If I need to see an example
# https://rdrr.io/github/chapmandu2/CollateralVulnerability2016/src/R/countHumanParalogs.R

