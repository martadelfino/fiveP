################################################################################
###### Script for creating counts of the gene paralogues #######################
################################################################################

library('tidyverse')



#' Calculate ratio for Paralogues
#'
#' @param paralogues A df of paralogue annotations
#' @param input_genes A vector of input genes
#' @return A dataframe with the fiveP scores
#' @export
calculate_paralogues_ratio <- function(paralogues, input_genes) {

  # Accessing the cleaned paralogues list from biomart ---------------------------

  #paralogues <- access_last_saved_file('./data/processed/', 'biomart_paralogues_cleaned') %>%
  #   dplyr::filter(complete.cases(paralog_hgnc_id))

#  rdata_paralogues <- load('data/biomart_paralogues_cleaned.RData')
 # paralogues <- get("biomart_paralogues_cleaned")

  rows_with_na <- paralogues %>%
    filter(is.na(paralog_hgnc_id))
 # print(rows_with_na)
  # there would have been 6603


  # Obtaining the input gene list ------------------------------------------------

  #input_genes <- data.frame(read_delim('./data/input/gene_list.txt', delim = '\t'))
  #print(input_genes)

#  rdata_input_genes <- load('data/input_genes.RData')
 # input_genes <- get("input_genes") %>%
  #  dplyr::select(hgnc_id)

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)


  # Remove anything below 30% - I've decided to take the mean --------------------

  paralogues_filtered <- paralogues %>%
    filter(mean_paralog_perc >= 30) %>%
    rename(hgnc_id = gene1_hgnc_id) %>%
    dplyr::select(hgnc_id, paralog_hgnc_id, mean_paralog_perc, max_paralog_perc)
 # head(paralogues_filtered)


  # Counting the number of paralogues per gene -----------------------------------

  paralogues_filtered_count1 <- paralogues_filtered %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(num_of_paralogs = n())
#  head(paralogues_filtered_count1)

  # Counting the number of paralogues that are input genes -----------------------

  # Check if paralog is an input gene
  paralogues_filtered_count2 <- paralogues_filtered_count1 %>%
    mutate(is_paralog_input_gene_yes_or_no = ifelse(paralog_hgnc_id %in% input_genes$hgnc_id, 1, 0))
#  print(paralogues_filtered_count2, n=5)

  # Count how many paralogues are input genes
  paralogues_filtered_count3 <- paralogues_filtered_count2 %>%
    group_by(hgnc_id) %>%
    mutate(num_input_gene_paralogs = sum(is_paralog_input_gene_yes_or_no))
#  print(paralogues_filtered_count3, n=5)

  # Save the current counts ------------------------------------------------------

  #save_data_with_datetime('./data/counts/', 'biomart_paralogues_intermediate_counts',
  #                       paralogues_filtered_count3)

#  save(paralogues_filtered_count3,
 #      file = "data/biomart_paralogues_intermediate_counts.RData")


  # Ratio of number of paralogs that are input gene : number of paralogs ---------

  paralogues_filtered_ratio <- paralogues_filtered_count3 %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(ratio_paraloginputgenes_to_paralogs = num_input_gene_paralogs / num_of_paralogs)
#  print(paralogues_filtered_ratio, n =5)

  paralogues_filtered_final <- paralogues_filtered_ratio %>%
    dplyr::select(!paralog_hgnc_id) %>%  dplyr::select(!mean_paralog_perc) %>%
    dplyr::select(!max_paralog_perc) %>%
    dplyr::distinct(hgnc_id, .keep_all = TRUE) %>%
    dplyr::select(!is_paralog_input_gene_yes_or_no) %>%
    dplyr::mutate(ratio_paraloginputgenes_to_paralogs = ifelse(is.na(ratio_paraloginputgenes_to_paralogs),
                                                               0, ratio_paraloginputgenes_to_paralogs))

#  print(paralogues_filtered_final, n=5)


  # Saving the results -----------------------------------------------------------

  #save_data_with_datetime('./data/counts', 'biomart_paralogues_ratio',
  #                       paralogues_filtered_final)

#  save(paralogues_filtered_final,
 #      file = "data/biomart_paralogues_ratio.RData")

  cat('\n(8/12) finished running paralogues_ratio.R\n')
  return(paralogues_filtered_final)

}
