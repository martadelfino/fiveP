################################################################################
###### Script for creating counts of the PPIs ##################################
################################################################################

library('tidyverse')


#' Calculate ratio for PPI
#'
#' @param ppi A df of ppi annotations
#' @param input_genes A df of input genes
#' @return A dataframe with the fiveP scores
#' @export
calculate_ppi_ratio <- function(ppi, input_genes) {
  # Accessing the PPIs from STRING -----------------------------------------------

  #ppi <- access_last_saved_file('./data/processed/', 'string_interactions') %>%
  # dplyr::select(!protein1_string_id) %>%
  #dplyr::select(!protein2_string_id) %>%
  #dplyr::rename(hgnc_id = protein1_hgnc_id)
  #head(ppi)

#  rdata_ppi <- load('data/string_interactions.RData')
 # ppi <- get("string_interactions") %>%
  #  dplyr::select(!protein1_string_id) %>%
   # dplyr::select(!protein2_string_id) %>%
    #dplyr::rename(hgnc_id = protein1_hgnc_id)

  ppi <- ppi %>%
    dplyr::select(!protein1_string_id) %>%
    dplyr::select(!protein2_string_id) %>%
    dplyr::rename(hgnc_id = protein1_hgnc_id)



  # Obtaining the input gene list ------------------------------------------------

  #input_genes <- data.frame(read_delim('./data/input/gene_list.txt', delim = '\t'))
  #print(input_genes)

#  rdata_input_genes <- load('data/input_genes.RData')
 # input_genes <- get("input_genes") %>%
  #  dplyr::select(hgnc_id)

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)

  # Counting the number of interactions per gene ---------------------------------

  ppi_count1 <- ppi %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(num_of_interactions = n())
 # head(ppi_count1)

  # Counting the number of interaction proteins that are input genes -------------

  # Check if interaction protein is an input gene
  ppi_count2 <- ppi_count1 %>%
    mutate(is_interaction_input_gene_yes_or_no = ifelse(protein2_hgnc_id %in% input_genes$hgnc_id, 1, 0))
#  print(ppi_count2, n=5)

  # Count how many paralogues are input genes
  ppi_count3 <- ppi_count2 %>%
    group_by(hgnc_id) %>%
    mutate(num_input_gene_interactions = sum(is_interaction_input_gene_yes_or_no))
#  print(ppi_count3, n=500)

  # Save the current counts ------------------------------------------------------

  #save_data_with_datetime('./data/counts/', 'string_interactions_intermediate_counts',
  #                       ppi_count3)

 # save(ppi_count3, file = "data/string_interactions_intermediate_counts.RData")


  # Ratio of number of interactors that are input gene : number of interactors ---

  ppi_ratio <- ppi_count3 %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(ratio_interactioninputgenes_to_interactions = num_input_gene_interactions / num_of_interactions)
#  print(ppi_ratio, n =5)

  ppi_ratio_final <- ppi_ratio %>%
    dplyr::select(!protein2_hgnc_id) %>%
    dplyr::distinct(hgnc_id, .keep_all = TRUE) %>%
    dplyr::select(!is_interaction_input_gene_yes_or_no) %>%
    dplyr::mutate(ratio_interactioninputgenes_to_interactions = ifelse(is.na(ratio_interactioninputgenes_to_interactions),
                                                                       0, ratio_interactioninputgenes_to_interactions))

#  print(ppi_ratio_final, n=50)


  # Saving the results -----------------------------------------------------------

  #save_data_with_datetime('./data/counts', 'string_interactions_ratio',
  #                       ppi_ratio_final)

 # save(ppi_ratio_final, file = "data/string_interactions_ratio.RData")

  print('finished running ppi_ratio.R')
  return(ppi_ratio_final)

}

