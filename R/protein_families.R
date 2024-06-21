################################################################################
###### Script for obtaining protein family participants ########################
################################################################################

library('tidyverse')

library('queryup')



#' Fetch Data
#'
#' This function fetches data from Database for a given list of genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param uniprot_input_gene_symbol_results_cleaned The df of uniprot results of input genes.
#' @return A dataframe with data from Database.
#' @export
fetch_protein_families <- function(protein_coding_genes,
                                   uniprot_input_gene_symbol_results_cleaned) {
  # Accessing the latest Uniprot results file ------------------------------------

  #uniprot_input_gene_symbol_results_cleaned <- access_last_saved_file('./data/processed/', 'uniprot_input_gene_symbol_results_cleaned')
  #print(uniprot_input_gene_symbol_results_cleaned$PANTHER)

#  rdata_uniprot_input_gene_symbol_results_cleaned <- load('data/uniprot_input_gene_symbol_results_cleaned.RData')
 # uniprot_input_gene_symbol_results_cleaned <- get("uniprot_input_gene_symbol_results_cleaned")

#  check_for_weird_characters(uniprot_input_gene_symbol_results_cleaned$PANTHER)


  # Creating a new df of the families --------------------------------------------

  input_genes_protein_families_expanded <- uniprot_input_gene_symbol_results_cleaned %>%
    tidyr::separate_rows(PANTHER, sep = ";") %>% distinct() %>%
    filter(PANTHER != "") %>%
    dplyr::select(PANTHER) %>% distinct() %>%
    dplyr::rename(family_id = PANTHER)
 # print(input_genes_protein_families_expanded, n=8)

  # removing extra information after the ':'
  input_genes_protein_families_expanded$family_id <- trimws(sub("\\:.*?", "", input_genes_protein_families_expanded$family_id))
 # print(input_genes_protein_families_expanded, n=10)


  # Querying -----------------------------------------------------------

  # Obtain families
  family <- dplyr::select(input_genes_protein_families_expanded, family_id)
#  print(family)

  # turning the object into a vector
  vector_family <- family %>% dplyr::pull(family_id)

  # Creating empty dataframe for the output
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
  uniprot_input_gene_family_results <- query_uniprot(query,
                                                     columns = columns,
                                                     show_progress = TRUE,
                                                     updateProgress = TRUE)

#  head(uniprot_input_gene_family_results)

  # Save the raw file ------------------------------------------------------------

  #save_data_with_datetime('./data/raw/', 'uniprot_pantherdb',
  #                       uniprot_input_gene_family_results)

  #save(uniprot_input_gene_family_results, file = "data/uniprot_pantherdb.RData")


  # Get the current date and time
  #current_datetime <- format(Sys.time(), "%Y_%m_%d__%H_%M_%S")

  # Create the file name with date and time
  #file_name <- paste0("./data/raw/uniprot_pantherdb_",
  #                   current_datetime, ".txt")

  # Save as a file
  #write.table(uniprot_input_gene_family_results, file_name, quote = F, sep = "\t", row.names = F)


  # Cleaning Protein Families result file from Uniprot ---------------------------

  # Selecting and renaming required columns
  proteinfamily_genes <- uniprot_input_gene_family_results %>%
    dplyr::select(Entry, HGNC, 'Gene Names (primary)', PANTHER) %>%
    dplyr::rename(uniprot_ids = Entry) %>%
    dplyr::rename(hgnc_id = HGNC) %>%
    dplyr::rename(family_id = PANTHER) %>%
    dplyr::rename(symbol = 'Gene Names (primary)')
 # print(head(proteinfamily_genes))

  # Removing trailing ;
  proteinfamily_genes$hgnc_id <- gsub(";$", "", proteinfamily_genes$hgnc_id)
  proteinfamily_genes$family_id <- gsub(";$", "", proteinfamily_genes$family_id)
#  head(proteinfamily_genes)

  # Separating families into new rows
  proteinfamily_genes_expanded <- proteinfamily_genes %>%
    tidyr::separate_rows(family_id, sep = ";") %>% # distinct() %>%
    filter(family_id != "")
#  print(proteinfamily_genes_expanded)

  # removing extra bits
  proteinfamily_genes_expanded$family_id <- trimws(sub("\\:.*", "", proteinfamily_genes_expanded$family_id))
  proteinfamily_genes_expanded <- proteinfamily_genes_expanded %>% distinct() %>%
    dplyr::select(family_id, uniprot_ids, hgnc_id, symbol) %>% # fixing order of columns
    arrange(family_id) # rearranging rows
#  print(proteinfamily_genes_expanded, n=100)



  # Saving the cleaned Protein Families file -------------------------------------

  #save_data_with_datetime('./data/processed/', 'uniprot_pantherdb_cleaned',
  #                       proteinfamily_genes_expanded)

  #save(proteinfamily_genes_expanded, file = "data/uniprot_pantherdb_cleaned.RData")

  print('finished running protein_families.R')
  return(proteinfamily_genes_expanded)

  # Calculating how many families are input gene families from total families-----

  # Access full protein families file from ftp
#  pantherdb_total_families <- read_delim('http://data.pantherdb.org/ftp/hmm_classifications/18.0/PANTHER18.0_HMM_classifications',
 #                                        delim = '\t', col_names = FALSE)
#  print(pantherdb_total_families)

  # Saving the raw file ----------------------------------------------------------

  #save_data_with_datetime('./data/raw/', 'pantherdb_total_families',
  #                       pantherdb_total_families)

  #save(pantherdb_total_families, file = "data/pantherdb_total_families.RData")


  # Cleaning the file ------------------------------------------------------------

  # Removing unnecessary columns and removing subfamilies

#  pantherdb_cleaned <- pantherdb_total_families %>%
 #   dplyr::select(X1) %>% rename(family_id = X1) %>%
 #   dplyr::mutate(family_id = str_replace(family_id, ":.*", "")) %>%
  #  distinct()
  # there are 15,963 families total

#  families_that_are_input_gene_families <- input_genes_protein_families_expanded %>%
 #   distinct()
  # 808 are NDD families



}

