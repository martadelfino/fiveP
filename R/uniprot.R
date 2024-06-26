################################################################################
###### Script to access Uniprot ################################################
################################################################################

library('queryup')

library('dplyr')

library('tidyverse')



#' Fetch Data
#'
#' This function fetches data from Database for a given list of genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param input_genes The df of input genes.
#' @return A dataframe with data from Database.
#' @export
fetch_uniprot <- function(protein_coding_genes, input_genes) {

  # Reading the protein coding genes file ----------------------------------------

  #hgnc_uniprot_symbol <- read_delim('./data/raw/protein_coding_genes.txt', '\t') %>%
  # dplyr::select(hgnc_id, uniprot_ids, symbol)
  #print(hgnc_uniprot_symbol)

 # rdata_hgnc_uniprot_symbol <- load('data/protein_coding_genes.RData')
#  hgnc_uniprot_symbol <- get("protein_coding_genes") %>%
  #  dplyr::select(hgnc_id, uniprot_ids, symbol)

  hgnc_uniprot_symbol <- protein_coding_genes %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol)


  # Checks of the protein coding genes file --------------------------------------

  # Checking for weird characters
#  check_for_weird_characters(hgnc_uniprot_symbol$hgnc_id)
 # check_for_weird_characters(hgnc_uniprot_symbol$uniprot_ids)
#  check_for_weird_characters(hgnc_uniprot_symbol$symbol)

  # checking NAs
  sum(is.na(hgnc_uniprot_symbol$uniprot_ids)) # 15 NAs
  sum(is.na(hgnc_uniprot_symbol$symbol)) # 0 NAs

  # Checking duplicates
  sum(duplicated(hgnc_uniprot_symbol$uniprot_ids)) # 128 duplicates

  # Removing all extra proteins identifiers (only keeping the canonical ones)
#  print(hgnc_uniprot_symbol$uniprot_ids)
  hgnc_uniprot_symbol$uniprot_ids <- trimws(sub('\\|.*', "", hgnc_uniprot_symbol$uniprot_ids))
#  print(hgnc_uniprot_symbol$uniprot_ids)
 # check_for_weird_characters(hgnc_uniprot_symbol$uniprot_ids)

  # Reading the input gene list --------------------------------------------------

  #input_genes <- read_delim('./data/input/gene_list.txt', delim = '\t')
  #print(input_genes)

#  rdata_input_genes <- load('data/input_genes.RData')
 # input_genes <- get("input_genes") %>%
  #  dplyr::select(hgnc_id)
  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)

  # Get the gene symbols for the input genes -------------------------------------

  input_genes_symbols <- input_genes %>% left_join(hgnc_uniprot_symbol,
                                                   by = 'hgnc_id')
#  print(input_genes_symbols)

  # Access uniprot ---------------------------------------------------------------

  # Note: Uniprot has no specific update schedule. So the entry version will be
  # important to keep track of database/entry versions.

  # Obtain gene symbols
  symbol <- dplyr::select(input_genes_symbols, symbol)
#  print(symbol)

  # turning the object into a vector
  vector_symbol <- symbol %>% dplyr::pull(symbol)

  # Creating empty dataframe for the output
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

  # Saving raw file --------------------------------------------------------------

  #save_data_with_datetime('./data/raw/', 'uniprot_input_gene_symbol_results_raw',
  #                       uniprot_input_gene_symbol_results)

  #save(uniprot_input_gene_symbol_results,
   #    file = "data/uniprot_input_gene_symbol_results_raw.RData")


  # Get the current date and time
  #current_datetime <- format(Sys.time(), "%Y_%m_%d__%H_%M_%S")

  # Create the file name with date and time
  #file_name <- paste0("./data/raw/uniprot_input_gene_symbol_results_raw_",
  #                   current_datetime, ".txt")

  # Save as a file
  #write.table(uniprot_input_gene_symbol_results, file_name,
  #            quote = F, sep = "\t", row.names = F)

  # Clean results ----------------------------------------------------------------

 # check_for_weird_characters(uniprot_input_gene_symbol_results)

  # Remove trailing ';' from the HGNC column
  uniprot_input_gene_symbol_results$HGNC <- gsub(";$", "",
                                                 uniprot_input_gene_symbol_results$HGNC)
 # head(uniprot_input_gene_symbol_results)

  # Expand any rows with multiple HGNCs
  rows_to_separate <- which(sapply(strsplit(uniprot_input_gene_symbol_results$HGNC, ";"),
                                   function(x) sum(grepl("HGNC:", x)) > 1))

  # View them
  cat('print(uniprot_input_gene_symbol_results[rows_to_separate,])')
  print(uniprot_input_gene_symbol_results[rows_to_separate,])

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
    #  head(uniprot_input_gene_symbol_results_combined)

  } else {
    uniprot_input_gene_symbol_results_combined <- uniprot_input_gene_symbol_results
  }


  # Checking duplicates
  sum(duplicated(uniprot_input_gene_symbol_results_combined$HGNC)) # 7 duplicates
  duplicated_rows <- uniprot_input_gene_symbol_results_combined[duplicated(uniprot_input_gene_symbol_results_combined), ]
#  print(duplicated_rows)

  # Remove duplicate rows
  uniprot_input_gene_symbol_results_combined <- distinct(uniprot_input_gene_symbol_results_combined)
#  head(uniprot_input_gene_symbol_results_combined)

  # checking duplicated HGNCs
  duplicated_rows_HGNC <- uniprot_input_gene_symbol_results_combined[duplicated(uniprot_input_gene_symbol_results_combined$HGNC), ]
#  print(duplicated_rows_HGNC)
  # some HGNCs just have more than one protein


  # Renaming columns
  uniprot_input_gene_symbol_results_cleaned <- uniprot_input_gene_symbol_results_combined %>%
    dplyr::rename(hgnc_id = HGNC) %>% dplyr::rename(uniprot_ids = Entry) %>%
    dplyr::rename(symbol = 'Gene Names (primary)')
#  head(uniprot_input_gene_symbol_results_cleaned)
#  print(duplicated(uniprot_input_gene_symbol_results_cleaned$uniprot_ids))



  merged_df <- merge(uniprot_input_gene_symbol_results_cleaned,
                     input_genes_symbols, by = c("hgnc_id", "uniprot_ids"))
#  head(merged_df)

  # removing extra columns after the merge
  uniprot_input_gene_symbol_results_cleaned <- merged_df %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol.x, ComplexPortal, PANTHER, 'Entry version') %>%
    dplyr::rename(symbol = symbol.x)
#  head(uniprot_input_gene_symbol_results_cleaned)

  #Specify the column to compare
  column_to_compare <- "hgnc_id"

  # Find rows where the values in the specified column are not the same
  not_same_rows <- anti_join(input_genes_symbols,
                             merged_df,
                             by = c(column_to_compare))

  # Print the result
 # print(not_same_rows)


  # Querying uniprot again in case there are missed proteins ---------------------

  ##### Create an if statement for this.


  # Obtain gene symbols
  symbol2 <- dplyr::select(not_same_rows, symbol)
 # print(symbol2)

  # turning the object into a vector
  vector_symbol2 <- symbol2 %>% dplyr::pull(symbol)

  # Creating empty dataframe for the output
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

#  print(uniprot_input_gene_symbol_results2)

  # Saving raw file of the missed proteins ---------------------------------------

  #save_data_with_datetime('./data/raw/', 'uniprot_input_gene_symbol_results_raw_missed_proteins',
  #                       uniprot_input_gene_symbol_results2)

  #save(uniprot_input_gene_symbol_results2,
   #    file = "data/uniprot_input_gene_symbol_results_raw_missed_proteins.RData")


  # Clean results ----------------------------------------------------------------

  # Remove trailing ';'
  uniprot_input_gene_symbol_results2$HGNC <- gsub(";$", "",
                                                  uniprot_input_gene_symbol_results2$HGNC)
#  head(uniprot_input_gene_symbol_results2)


  ######### create another if statement to check if there are any rows to separate ------------
  # Expand any rows with multiple HGNCs
  rows_to_separate2 <- which(sapply(strsplit(uniprot_input_gene_symbol_results2$HGNC, ";"),
                                    function(x) sum(grepl("HGNC:", x)) > 1))

  if (length(rows_to_separate2) > 0) {
    # View them
 #   print(uniprot_input_gene_symbol_results2[rows_to_separate2,])

    uniprot_input_gene_symbol_results_separated2 <- uniprot_input_gene_symbol_results2 %>%
      slice(rows_to_separate2) %>%
      separate_rows(HGNC, sep = ";") %>%
      mutate(HGNC = ifelse(grepl("HGNC:", HGNC), HGNC, paste0("HGNC:", HGNC))) %>%
      as.data.frame() %>%
      mutate_if(is.logical, as.character)
 #   print(uniprot_input_gene_symbol_results_separated2)

    # Combine the separated rows with the rest of the df
    uniprot_input_gene_symbol_results_combined2 <- bind_rows(uniprot_input_gene_symbol_results2[-rows_to_separate2, ],
                                                             uniprot_input_gene_symbol_results_separated2)
 #   head(uniprot_input_gene_symbol_results_combined2)

    # Checking duplicates
    sum(duplicated(uniprot_input_gene_symbol_results_combined2$HGNC)) # 0 duplicates
    duplicated_rows2 <- uniprot_input_gene_symbol_results_combined2[duplicated(uniprot_input_gene_symbol_results_combined2), ]
 #   print(duplicated_rows2)

    # Remove duplicate rows
    uniprot_input_gene_symbol_results_combined2 <- distinct(uniprot_input_gene_symbol_results_combined2)
#    head(uniprot_input_gene_symbol_results_combined2)

    # checking duplicated HGNCs
    duplicated_rows_HGNC2 <- uniprot_input_gene_symbol_results_combined[duplicated(uniprot_input_gene_symbol_results_combined2$HGNC), ]
#    print(duplicated_rows_HGNC2)

  } else {
    print('No duplicated rows found. Continuing to the next step.')
  }


  # If there are no rows to separate, continue straight to renaming columns -------------------

  # Renaming columns
  uniprot_input_gene_symbol_results2_cleaned <- uniprot_input_gene_symbol_results2 %>%
    dplyr::rename(hgnc_id = HGNC) %>% dplyr::rename(uniprot_ids = Entry) %>%
    dplyr::rename(symbol = 'Gene Names (primary)')
#  head(uniprot_input_gene_symbol_results2_cleaned)


  # Adding these rows to the first results ---------------------------------------

  final_uniprot_input_genes_results <- bind_rows(uniprot_input_gene_symbol_results_cleaned,
                                                 uniprot_input_gene_symbol_results2_cleaned)
#  tail(final_uniprot_input_genes_results)


  # Save the FINALLLLLLL cleaned file --------------------------------------------------------

  #save_data_with_datetime('./data/processed/', 'uniprot_input_gene_symbol_results_cleaned',
  #                       final_uniprot_input_genes_results)

  #save(final_uniprot_input_genes_results,
   #    file = "data/uniprot_input_gene_symbol_results_cleaned.RData")

  cat('\n(5/12) finished running uniprot.R\n')
  return(final_uniprot_input_genes_results)
}


