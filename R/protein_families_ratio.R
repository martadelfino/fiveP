################################################################################
###### Script for creating counts of the protein families ######################
################################################################################

library('tidyverse')

#source('./src/functions.R')



calculate_protein_families_ratio <- function(panther, input_genes) {
  # Accessing the cleaned gene list queried Uniprot PantherDB file ---------------

  #panther <- access_last_saved_file('./data/processed/', 'uniprot_pantherdb_cleaned')

 # rdata_panther <- load('data/uniprot_pantherdb_cleaned.RData')
#  panther <- get("uniprot_pantherdb_cleaned")

 # check_for_weird_characters(panther$family_id)

  # Counting number of input proteins/genes in each family -----------------------

  #input_genes <- read_delim('./data/input/gene_list.txt', delim = '\t')
  #print(input_genes)

#  rdata_input_genes <- load('data/input_genes.RData')
 # input_genes <- get("input_genes") %>%
  #  dplyr::select(hgnc_id)

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)


  panther_counts <- panther %>%
    mutate(input_gene_yes_or_no = ifelse(hgnc_id %in% input_genes$hgnc_id, 1, 0))
 # print(panther_counts, n=5)

  # Counting the number of input genes per family

  panther_counts <- panther_counts %>%
    group_by(family_id) %>%
    mutate(num_genes_in_family = n(),
           num_input_gene_per_family = sum(input_gene_yes_or_no))
 # print(panther_counts, n=5)

  # Saving the counts ------------------------------------------------------------

  #save_data_with_datetime('./data/counts/', 'uniprot_pantherdb_counts',
  #                       panther_counts)

 # save(panther_counts, file = "data/uniprot_pantherdb_counts.RData")



  # Creating a df for the number of unique genes in each family ------------------

  # Counting the number of unique genes in each pathway that gene is related to
  panther_counts_per_gene <- panther_counts %>%
    group_by(hgnc_id) %>%
    dplyr::mutate(
      num_families = n_distinct(family_id),
      num_unique_genes_in_families = sum(length(unique(panther_counts$hgnc_id[panther_counts$family_id %in% family_id])) - 1),
      num_input_genes_in_families = sum(unique(panther_counts$hgnc_id[panther_counts$family_id %in% family_id]) %in% input_genes$hgnc_id) - (hgnc_id %in% input_genes$hgnc_id)
    ) %>%
    dplyr::mutate(ratio_input_genes_in_families = num_input_genes_in_families / num_unique_genes_in_families) %>%
    dplyr::select(hgnc_id, uniprot_ids, family_id, num_families,
                  num_unique_genes_in_families, num_input_genes_in_families,
                  ratio_input_genes_in_families) %>% arrange(hgnc_id)
 # print(panther_counts_per_gene, n= 10)

  panther_counts_per_gene_final <- panther_counts_per_gene %>%
    dplyr::select(hgnc_id, uniprot_ids, num_families,
                  num_unique_genes_in_families, num_input_genes_in_families,
                  ratio_input_genes_in_families) %>% unique() %>%
    dplyr::mutate(ratio_input_genes_in_families = ifelse(is.na(ratio_input_genes_in_families), 0, ratio_input_genes_in_families))
 # print(panther_counts_per_gene_final, n= 10)

  # Saving the pathway counts per gene -------------------------------------------

  #save_data_with_datetime('./data/counts/', 'uniprot_pantherdb_gene_counts',
  #                       panther_counts_per_gene_final)

#  save(panther_counts_per_gene_final, file = "data/uniprot_pantherdb_gene_counts.RData")

  print('finished running protein_families_ratio.R')
  return(panther_counts_per_gene_final)


}
