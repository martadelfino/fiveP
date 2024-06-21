################################################################################
###### Script for obtaining protein protein interactions #######################
################################################################################

library("tidyverse")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("STRINGdb")
library('STRINGdb')
packageVersion('STRINGdb')

#' Fetch Data
#'
#' This function fetches data from Database for a given list of genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @return A dataframe with data from Database.
#' @export

fetch_ppi <- function(protein_coding_genes) {

  # Reading the protein coding genes file ----------------------------------------

  #hgnc_esembl <- read_delim('./data/raw/protein_coding_genes.txt', '\t') %>%
  # dplyr::select(hgnc_id, ensembl_gene_id)
  #print(hgnc_esembl)

#  rdata_hgnc_esembl <- load('data/protein_coding_genes.RData')
 # hgnc_esembl <- get("protein_coding_genes") %>%
  #  dplyr::select(hgnc_id, ensembl_gene_id)

  hgnc_esembl <- protein_coding_genes %>%
    dplyr::select(hgnc_id, ensembl_gene_id)

  hgnc_ensembl <- data.frame(hgnc_esembl)

  # Load STRING database ---------------------------------------------------------

  #string_db <- STRINGdb$new(version = "12.0", species = 9606,
  #                        score_threshold = 700, input_directory="")

  string_db <- STRINGdb$new(version="12.0", species=9606,
                            score_threshold=700, network_type="full",
                            input_directory="")
  class(string_db)

  # Mapping input gene list to the STRING identifiers ----------------------------

  input_genes_mapped <- string_db$map(hgnc_ensembl, "ensembl_gene_id", removeUnmappedRows = TRUE )
  # Warning:  we couldn't map to STRING 0% of your identifiers
 # print(input_genes_mapped)

  input_genes_mapped_vector <- input_genes_mapped %>% pull(STRING_id)


  # Get interactions -------------------------------------------------------------

  interactions <- string_db$get_interactions(input_genes_mapped_vector)
 # head(interactions)

  min(interactions$combined_score)

  # remove duplicate ones
  interactions <- interactions %>%
    dplyr::distinct()

  # Saving the STRING PPI interactions file ----------------------------------

  #save_data_with_datetime('./data/raw', 'string_interactions', interactions)

  #save(interactions, file = "data/string_interactions_raw.RData")


  # Clean the file ---------------------------------------------------------------

  # Map back to HGNC IDs
  interactions_cleaned <- input_genes_mapped %>%
    left_join(interactions, join_by(STRING_id == from),
              relationship = "many-to-many") %>%
    dplyr::select(hgnc_id, STRING_id, to, combined_score) %>%
    dplyr::rename(protein1_hgnc_id = hgnc_id) %>%
    dplyr::rename(protein1_string_id = STRING_id) %>%
    dplyr::rename(protein2_string_id = to)
#  head(interactions_cleaned)

  interactions_cleaned_hgnc <- input_genes_mapped %>%
    left_join(interactions_cleaned, join_by(STRING_id == protein2_string_id),
              relationship = "many-to-many") %>%
    dplyr::select(protein1_hgnc_id, protein1_string_id, hgnc_id, STRING_id, combined_score) %>%
    dplyr::rename(protein2_hgnc_id = hgnc_id) %>%
    dplyr::rename(protein2_string_id = STRING_id) %>%
    dplyr::arrange(protein1_hgnc_id)

#  head(interactions_cleaned_hgnc)

  # Saving the processed STRING PPI interactions file ----------------------------

  #save_data_with_datetime('./data/processed', 'string_interactions',
  #                       interactions_cleaned_hgnc)

  #save(interactions_cleaned_hgnc, file = "data/string_interactions.RData")

  print('finished running ppi.R')
  return(interactions_cleaned_hgnc)

}
