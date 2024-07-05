library('queryup')

library('dplyr')

library('tidyverse')


#' Fetch Uniprot Data
#'
#' This function fetches protein data from Uniprot for all protein coding genes_______?.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param input_genes The df of input genes.
#' @return A dataframe with Uniprot data of the input genes.
#' @export
fetch_uniprot <- function(protein_coding_genes, input_genes) {

  # Reading the protein coding genes file --------------------------------------

  hgnc_uniprot_symbol <- protein_coding_genes %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol)


  # Checks of the protein coding genes file ------------------------------------

  # Removing all extra proteins identifiers (only keeping the canonical ones)
  hgnc_uniprot_symbol$uniprot_ids <- trimws(sub('\\|.*', "", hgnc_uniprot_symbol$uniprot_ids))


  # Reading the input gene list --------------------------------------------------

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)

  # Get gene symbols for the input genes
  input_genes_symbols <- input_genes %>% left_join(hgnc_uniprot_symbol,
                                                   by = 'hgnc_id')

  # Access uniprot ---------------------------------------------------------------

  # Note: Uniprot has no specific update schedule. So the entry version will be
  # important to keep track of database/entry versions.

  # Obtain gene symbols
  symbol <- dplyr::select(input_genes_symbols, symbol)
  vector_symbol <- symbol %>% dplyr::pull(symbol)  # turn object into vector

  # Create empty df for the output
  uniprot_input_gene_symbol_results <- data.frame()

  # Construct the query
  query <- list("gene_exact" = vector_symbol,
                "reviewed" = "true",
                "organism_id" = "9606")
  columns <- c('id',
               'xref_hgnc',
               'gene_primary',
               "xref_complexportal",
               'xref_panther',
               'version')
  uniprot_input_gene_symbol_results <- queryup::query_uniprot(query,
                                                              columns = columns,
                                                              show_progress = TRUE,
                                                              updateProgress = TRUE)


  ## Clean results
  # Remove trailing ';' from the HGNC column
  uniprot_input_gene_symbol_results$HGNC <- gsub(";$", "",
                                                 uniprot_input_gene_symbol_results$HGNC)
  # Expand any rows with multiple HGNCs
  rows_to_separate <- which(sapply(strsplit(uniprot_input_gene_symbol_results$HGNC, ";"),
                                   function(x) sum(grepl("HGNC:", x)) > 1))
  # Check if rows_to_separate is empty
  if (length(rows_to_separate) > 0) {
    uniprot_input_gene_symbol_results_separated <- uniprot_input_gene_symbol_results %>%
      slice(rows_to_separate) %>%
      separate_rows(HGNC, sep = ";") %>%
      mutate(HGNC = ifelse(grepl("HGNC:", HGNC), HGNC, paste0("HGNC:", HGNC)))
    cat('print(uniprot_input_gene_symbol_results_separated)')
    print(uniprot_input_gene_symbol_results_separated)

    # Combine the separated rows with the rest of the dataframe
    uniprot_input_gene_symbol_results_combined <- bind_rows(uniprot_input_gene_symbol_results[-rows_to_separate, ],
                                                            uniprot_input_gene_symbol_results_separated)

  } else { # continue
    uniprot_input_gene_symbol_results_combined <- uniprot_input_gene_symbol_results
  }

  # Checking duplicates
  sum(duplicated(uniprot_input_gene_symbol_results_combined$HGNC))
  duplicated_rows <- uniprot_input_gene_symbol_results_combined[duplicated(uniprot_input_gene_symbol_results_combined), ]
  # Remove duplicate rows
  uniprot_input_gene_symbol_results_combined <- distinct(uniprot_input_gene_symbol_results_combined)

  # checking duplicated HGNCs
  duplicated_rows_HGNC <- uniprot_input_gene_symbol_results_combined[duplicated(uniprot_input_gene_symbol_results_combined$HGNC), ]
  # some HGNCs just have more than one protein

  # Renaming columns
  uniprot_input_gene_symbol_results_cleaned <- uniprot_input_gene_symbol_results_combined %>%
    dplyr::rename(hgnc_id = HGNC) %>% dplyr::rename(uniprot_ids = Entry) %>%
    dplyr::rename(symbol = 'Gene Names (primary)')


  # Merging the cleaned data with the input gene symbols -----------------------

  merged_df <- merge(uniprot_input_gene_symbol_results_cleaned,
                     input_genes_symbols, by = c("hgnc_id", "uniprot_ids"))

  # removing extra columns after the merge
  uniprot_input_gene_symbol_results_cleaned <- merged_df %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol.x, ComplexPortal, PANTHER, 'Entry version') %>%
    dplyr::rename(symbol = symbol.x)


  # Querying uniprot again in case there are missed proteins ---------------------

  #Specify the column to compare
  column_to_compare <- "hgnc_id"

  # Find rows where the values in the specified column are not the same
  not_same_rows <- anti_join(input_genes_symbols,
                             merged_df,
                             by = c(column_to_compare))

  ##### Create an if statement for this.

  # Obtain gene symbols
  symbol2 <- dplyr::select(not_same_rows, symbol)
  vector_symbol2 <- symbol2 %>% dplyr::pull(symbol)
  uniprot_input_gene_symbol_results2 <- data.frame()

  # Construct the query
  query2 <- list("gene_exact" = vector_symbol2,
                 "reviewed" = "true",
                 "organism_id" = "9606")
  columns2 <- c('id',
                'xref_hgnc',
                'gene_primary',
                "xref_complexportal",
                'xref_panther',
                'version')
  uniprot_input_gene_symbol_results2 <- queryup::query_uniprot(query2,
                                                               columns = columns,
                                                               show_progress = TRUE,
                                                               updateProgress = TRUE)


  # Clean the results

  # Remove trailing ';'
  uniprot_input_gene_symbol_results2$HGNC <- gsub(";$", "",
                                                  uniprot_input_gene_symbol_results2$HGNC)

  # Expand any rows with multiple HGNCs
  rows_to_separate2 <- which(sapply(strsplit(uniprot_input_gene_symbol_results2$HGNC, ";"),
                                    function(x) sum(grepl("HGNC:", x)) > 1))

  if (length(rows_to_separate2) > 0) {

    uniprot_input_gene_symbol_results_separated2 <- uniprot_input_gene_symbol_results2 %>%
      slice(rows_to_separate2) %>%
      separate_rows(HGNC, sep = ";") %>%
      mutate(HGNC = ifelse(grepl("HGNC:", HGNC), HGNC, paste0("HGNC:", HGNC))) %>%
      as.data.frame() %>%
      mutate_if(is.logical, as.character)

    # Combine the separated rows with the rest of the df
    uniprot_input_gene_symbol_results_combined2 <- bind_rows(uniprot_input_gene_symbol_results2[-rows_to_separate2, ],
                                                             uniprot_input_gene_symbol_results_separated2)
    # Checking duplicates
    sum(duplicated(uniprot_input_gene_symbol_results_combined2$HGNC)) # 0 duplicates
    duplicated_rows2 <- uniprot_input_gene_symbol_results_combined2[duplicated(uniprot_input_gene_symbol_results_combined2), ]

    # Remove duplicate rows
    uniprot_input_gene_symbol_results_combined2 <- distinct(uniprot_input_gene_symbol_results_combined2)

    # checking duplicated HGNCs
    duplicated_rows_HGNC2 <- uniprot_input_gene_symbol_results_combined[duplicated(uniprot_input_gene_symbol_results_combined2$HGNC), ]

  } else {
    print('No duplicated rows found. Continuing to the next step.')
  }


  # If there are no rows to separate, continue straight to renaming columns -------------------

  # Renaming columns
  uniprot_input_gene_symbol_results2_cleaned <- uniprot_input_gene_symbol_results2 %>%
    dplyr::rename(hgnc_id = HGNC) %>% dplyr::rename(uniprot_ids = Entry) %>%
    dplyr::rename(symbol = 'Gene Names (primary)')

  # Adding these rows to the first results

  final_uniprot_input_genes_results <- bind_rows(uniprot_input_gene_symbol_results_cleaned,
                                                 uniprot_input_gene_symbol_results2_cleaned)


  cat('\n(5/12) finished running uniprot.R\n')
  return(final_uniprot_input_genes_results)
}


