library('tidyverse')

library('queryup')


#' Fetch Protein Family Data from Uniprot
#'
#' This function fetches protein family data from Uniprot for the input genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param uniprot_input_gene_symbol_results_cleaned The df of uniprot results of input genes.
#' @return A dataframe with Protein Family data from Uniprot for the input genes.
#' @export
fetch_protein_families <- function(protein_coding_genes,
                                   uniprot_input_gene_symbol_results_cleaned) {

  # Creating a df of protein families data from uniprot results ----------------

  input_genes_protein_families_expanded <- uniprot_input_gene_symbol_results_cleaned %>%
    tidyr::separate_rows(PANTHER, sep = ";") %>% distinct() %>%
    filter(PANTHER != "") %>%
    dplyr::select(PANTHER) %>% distinct() %>%
    dplyr::rename(family_id = PANTHER)

  # removing extra information after the ':'
  input_genes_protein_families_expanded$family_id <- trimws(sub("\\:.*?", "", input_genes_protein_families_expanded$family_id))


  # Querying Uniprot -----------------------------------------------------------

  # Obtain families
  family <- dplyr::select(input_genes_protein_families_expanded, family_id)
  vector_family <- family %>% dplyr::pull(family_id)   # turning object into vector

  # Empty dataframe for the output
  uniprot_input_gene_family_results <- data.frame()

  # Construct the query
  query <- list("xref" = vector_family,
                "reviewed" = "true",
                "organism_id" = "9606")
  columns <- c('id',
               'xref_hgnc',
               'gene_primary',
               'xref_panther',
               'version')
  uniprot_input_gene_family_results <- queryup::query_uniprot(query,
                                                              columns = columns,
                                                              show_progress = TRUE,
                                                              updateProgress = TRUE)


  # Cleaning Protein Families result file from Uniprot ---------------------------

  # Selecting and renaming required columns
  proteinfamily_genes <- uniprot_input_gene_family_results %>%
    dplyr::select(Entry, HGNC, 'Gene Names (primary)', PANTHER) %>%
    dplyr::rename(uniprot_ids = Entry) %>%
    dplyr::rename(hgnc_id = HGNC) %>%
    dplyr::rename(family_id = PANTHER) %>%
    dplyr::rename(symbol = 'Gene Names (primary)')

  # Removing trailing ;
  proteinfamily_genes$hgnc_id <- gsub(";$", "", proteinfamily_genes$hgnc_id)
  proteinfamily_genes$family_id <- gsub(";$", "", proteinfamily_genes$family_id)

  # Separating families into new rows
  proteinfamily_genes_expanded <- proteinfamily_genes %>%
    tidyr::separate_rows(family_id, sep = ";") %>%
    filter(family_id != "")

  # removing extra bits
  proteinfamily_genes_expanded$family_id <- trimws(sub("\\:.*", "", proteinfamily_genes_expanded$family_id))
  proteinfamily_genes_expanded <- proteinfamily_genes_expanded %>% distinct() %>%
    dplyr::select(family_id, uniprot_ids, hgnc_id, symbol) %>% # fixing order of columns
    arrange(family_id) # rearranging rows


  cat('\n(7/12) finished running protein_families.R\n')
  return(proteinfamily_genes_expanded)

}

