################################################################################
###### Script for creating counts of the protein complexes #####################
################################################################################

library('tidyverse')


#' Calculate ratio for Protein Complex
#'
#' @param complexportal A df of protein complex annotations
#' @param input_genes A vector of input genes
#' @return A dataframe with the fiveP scores
#' @export
calculate_protein_complex_ratio <- function(complexportal, input_genes) {
  # Accessing the cleaned gene list queried ComplexPortal file -------------------

  #complexportal <- access_last_saved_file('./data/processed/', 'complexportal_cleaned')

 # rdata_complexportal <- load('data/complexportal_cleaned.RData')
  #complexportal <- get("complexportal_cleaned")

#  check_for_weird_characters(complexportal$complex_id)


  # Counting number of input proteins/genes in each complex ----------------------

  #input_genes <- read_delim('./data/input/gene_list.txt', delim = '\t')
  #print(input_genes)

#  rdata_input_genes <- load('data/input_genes.RData')
 # input_genes <- get("input_genes") %>%
  #  dplyr::select(hgnc_id)

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)

  complexportal_counts <- complexportal %>%
    mutate(input_gene_yes_or_no = ifelse(hgnc_id %in% input_genes$hgnc_id, 1, 0))
#  print(complexportal_counts, n=5)

  # Counting the number of input genes per complex

  complexportal_counts <- complexportal_counts %>%
    group_by(complex_id) %>%
    mutate(num_genes_in_complex = n(),
           num_input_gene_per_complex = sum(input_gene_yes_or_no))
#  print(complexportal_counts, n=5)

  # Saving the counts ------------------------------------------------------------

  #save_data_with_datetime('./data/counts/', 'complexportal_counts',
  #                       complexportal_counts)

#  save(complexportal_counts, file = "data/complexportal_counts.RData")


  # Creating a df for the number of unique genes in each protein complex ---------

  # Counting the number of unique genes in each complex that gene is related to
  complexportal_counts_per_gene <- complexportal_counts %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(
      num_complexes = n_distinct(complex_id),
      num_unique_genes_in_complexes = sum(length(unique(complexportal_counts$hgnc_id[complexportal_counts$complex_id %in% complex_id])) - 1),
      num_input_genes_in_complexes = sum(unique(complexportal_counts$hgnc_id[complexportal_counts$complex_id %in% complex_id]) %in% input_genes$hgnc_id) - (hgnc_id %in% input_genes$hgnc_id)
    ) %>%
    dplyr::mutate(ratio_input_genes_in_complexes = num_input_genes_in_complexes / num_unique_genes_in_complexes) %>%
    dplyr::select(hgnc_id, symbol, uniprot_ids, complex_id, num_complexes,
                  num_unique_genes_in_complexes, num_input_genes_in_complexes,
                  ratio_input_genes_in_complexes)
#  print(complexportal_counts_per_gene, n= 100)

  complexportal_counts_per_gene_final <- complexportal_counts_per_gene %>%
    dplyr::select(hgnc_id, symbol, uniprot_ids, num_complexes,
                  num_unique_genes_in_complexes, num_input_genes_in_complexes,
                  ratio_input_genes_in_complexes) %>% unique() %>%
    dplyr::mutate(ratio_input_genes_in_complexes = ifelse(is.na(ratio_input_genes_in_complexes),
                                                          0, ratio_input_genes_in_complexes))
#  print(complexportal_counts_per_gene_final, n= 10)


  # Saving the pathway counts per gene -------------------------------------------

  #save_data_with_datetime('./data/counts/', 'complexportal_gene_counts',
  #                       complexportal_counts_per_gene_final)

 # save(complexportal_counts_per_gene_final, file = "data/complexportal_gene_counts.RData")

  cat('\n(11/12) finished running protein_complex_ratio.R\n')
  return(complexportal_counts_per_gene_final)

}
