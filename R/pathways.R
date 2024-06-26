################################################################################
###### Script for obtaining pathways ###########################################
################################################################################

library('tidyverse')



#' Fetch Data
#'
#' This function fetches data from Database for a given list of genes.
#'
#' @param protein_coding_genes The df of all protein coding genes.
#' @param input_genes The df of input genes.
#' @return A list of two dataframes with data from Database.
#' @export
fetch_pathways <- function(protein_coding_genes, input_genes) {
  # Reading the protein coding genes file ----------------------------------------

  #hgnc_uniprot_symbol_entrez <- read_delim('./data/raw/protein_coding_genes.txt', '\t') %>%
  # dplyr::select(hgnc_id, uniprot_ids, symbol, entrez_id)
  #print(hgnc_uniprot_symbol_entrez)

#  rdata_hgnc_uniprot_symbol_entrez <- load('data/protein_coding_genes.RData')
 # hgnc_uniprot_symbol_entrez <- get("protein_coding_genes") %>%
  #  dplyr::select(hgnc_id, uniprot_ids, symbol, entrez_id)

  hgnc_uniprot_symbol_entrez <- protein_coding_genes %>%
    dplyr::select(hgnc_id, uniprot_ids, symbol, entrez_id)

  # Checks of the protein coding genes file
 # check_for_weird_characters(hgnc_uniprot_symbol_entrez$entrez_id)

  # Removing all extra proteins identifiers (only keeping the canonical ones)
  #print(hgnc_uniprot_symbol_entrez$uniprot_ids)
  hgnc_uniprot_symbol_entrez$uniprot_ids <- trimws(sub('\\|.*', "", hgnc_uniprot_symbol_entrez$uniprot_ids))
  #print(hgnc_uniprot_symbol_entrez$uniprot_ids)


  # Reading the input gene list --------------------------------------------------

  #input_genes <- read_delim('./data/input/gene_list.txt', delim = '\t')
  #print(input_genes)

#  rdata_input_genes <- load('data/input_genes.RData')
 # input_genes <- get("input_genes") %>%
  #  dplyr::select(hgnc_id)

  input_genes <- input_genes %>%
    dplyr::select(hgnc_id)


  # Get the gene symbols for the input genes -------------------------------------

  input_genes_entrez_id <- input_genes %>%
    inner_join(hgnc_uniprot_symbol_entrez, by = 'hgnc_id') %>%
    pull(entrez_id) %>%
    as.character(.)

 # print(input_genes_entrez_id)


  # Reactome enriched pathways ---------------------------------------------------

  #input_genes_reactome_enriched <- enrichPathway(gene = input_genes_entrez_id,
  #                                      pvalueCutoff = 0.05,
  #                                      readable = TRUE)
  #head(input_genes_reactome_enriched)

  # Get the summary of enriched pathways per gene
  #input_genes_reactome_enriched_df <- as.data.frame(input_genes_reactome_enriched)
  #head(input_genes_reactome_enriched_df)

  #unique_test <- data.frame(input_genes_reactome_enriched_df) %>%
  #  count(ID)
  #print(unique_test)


  # Saving the raw enriched Reactome pathways file -------------------------------

  #save_data_with_datetime("./data/raw/", "reactome_enriched",
  #                        input_genes_reactome_enriched_df)


  # Save the plot of the results -------------------------------------------------

  #input_genes_reactome_pairwise_sim <- pairwise_termsim(input_genes_reactome_enriched)

  # Get the current date and time
  #current_datetime_plot <- format(Sys.time(), "%Y_%m_%d__%H_%M_%S")
  # Create the file name with date and time
  #file_name_plot <- paste0("./data/raw/plot_reactome_enriched_",
  #                   current_datetime_plot, ".png")
  #png(file_name_plot, width = 6000, height = 4000, res = 300)
  #emapplot(input_genes_reactome_pairwise_sim)
  #dev.off()



  # Retrieving the genes in all pathways -----------------------------------------

  # Obtaining data from Reactome directly, lowest level pathways
  Uniprot2Reactome <- read_delim('https://reactome.org/download/current/UniProt2Reactome.txt',
                                 col_names = FALSE)
 # head(Uniprot2Reactome)

  # Saving the raw Uniprot to Reactome pathway file ------------------------------

  #save_data_with_datetime("./data/raw/", "reactome_Uniprot2Reactome",
  #                       Uniprot2Reactome)

  #save(Uniprot2Reactome, file = "data/reactome_Uniprot2Reactome.RData")



  # Cleaning the Uniprot to Reactome file ----------------------------------------

  # renaming columns
  Uniprot2Reactome_cleaned <- Uniprot2Reactome %>%
    dplyr::rename(uniprot_ids = X1, pathway_id = X2, url = X3, pathway_name = X4,
                  evidence_code = X5, species = X6)
 # head(Uniprot2Reactome_cleaned)

  # selecting homo sapiens only
  Uniprot2Reactome_cleaned <- Uniprot2Reactome_cleaned %>%
    dplyr::filter(species == 'Homo sapiens')
#  head(Uniprot2Reactome_cleaned)

  # selecting Uniprot id and Pathway id columns only
  Uniprot2Reactome_final <- Uniprot2Reactome_cleaned %>%
    dplyr::select(uniprot_ids, pathway_id)
 # head(Uniprot2Reactome_final)

  # Joining the file with the protein coding genes file
  Uniprot2Reactome_final_hgnc <- Uniprot2Reactome_final %>%
    left_join(hgnc_uniprot_symbol_entrez, by = 'uniprot_ids', relationship = 'many-to-many')
#  print(Uniprot2Reactome_final_hgnc, n = 100)
  # Some proteins don't map to hgnc ids. this is because they are immunoglobulins
  # or other immune related proteins. These are ignore.

  Uniprot2Reactome_final_hgnc_no_na <- Uniprot2Reactome_final_hgnc %>%
    filter(!is.na(hgnc_id))
 # print(Uniprot2Reactome_final_hgnc_no_na, n=100)

  unique_test <- data.frame(Uniprot2Reactome_final_hgnc_no_na) %>%
    count(hgnc_id)
#  print(unique_test) # 10959 genes

  # Saving the cleaned Uniprot to Reactome pathway file --------------------------

  #save_data_with_datetime("./data/processed/", "reactome_Uniprot2Reactome_cleaned",
  #                       Uniprot2Reactome_final_hgnc_no_na)

  #save(Uniprot2Reactome_final_hgnc_no_na, file = "data/reactome_Uniprot2Reactome_cleaned.RData")


  # Creating a file with only the pathways of the input genes --------------------

  # Obtain uniprot ids for the input genes
  #input_genes_uniprot_ids <- input_genes %>%
  #  left_join(hgnc_uniprot_symbol, by = 'hgnc_id') %>%
  #  dplyr::select(uniprot_ids)
  #print(input_genes_uniprot_ids)

  # Filter the Uniprot2Reactome final file with the input genes
  input_genes_Uniprot2Reactome <- Uniprot2Reactome_final_hgnc_no_na %>%
    dplyr::filter(hgnc_id %in% input_genes$hgnc_id)
#  print(input_genes_Uniprot2Reactome)

  unique_test <- data.frame(input_genes_Uniprot2Reactome) %>%
    count(pathway_id)
#  print(unique_test)


  # Saving the input gene list pathways Uniprot to Reactome file -----------------

  #save_data_with_datetime("./data/processed/", "reactome_gene_list_Uniprot2Reactome",
  #                       input_genes_Uniprot2Reactome)

  #save(input_genes_Uniprot2Reactome, file = "data/reactome_gene_list_Uniprot2Reactome.RData")

  cat('\n(3/12) finished running pathways.R\n')
  return(list(input_genes_Uniprot2Reactome = input_genes_Uniprot2Reactome,
              Uniprot2Reactome_final_hgnc_no_na = Uniprot2Reactome_final_hgnc_no_na))

  # Calculating how many pathways are input gene list pathways from the total pathways ----

  #reactome_total_pathways <- Uniprot2Reactome_final %>%
   # dplyr::select(pathway_id) %>% distinct()
  # there are 2167 pathways in total

  #pathways_that_are_input_gene_pathways <- input_genes_Uniprot2Reactome %>%
   # dplyr::select(pathway_id) %>% distinct()
  #  878 are NDD pathways
}



