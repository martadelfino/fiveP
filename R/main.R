# main file. still to do


# Read the input gene list. Choose which column in dataframe I need.

# must be a df with a column called 'hgnc_id'
read_input_genes <- function(df) {
  df <- df %>%
    dplyr::select(hgnc_id)
  save(df, file = "data/input_genes.RData")
}


# Run the data fetching scripts

fetch_data <- function() {
  script_paths <- c(
    "R/hgnc_gene_list.R",
    "R/uniprot.R",
    "R/paralogues.R",
    "R/pathways.R",
    "R/ppi.R",
    "R/protein_complex.R",
    "R/protein_families.R"
  )

  for (script_path in script_paths) {
    source(script_path)
  }
}


#test <- fetch_data()

# Run the count calculations

create_scores <- function() {
  script_paths <- c(
    "R/paralogues_ratio.R",
    "R/pathways_ratio.R",
    "R/ppi_ratio.R",
    "R/protein_complex_ratio.R",
    "R/protein_families_ratio.R"
  )

  for (script_path in script_paths) {
    source(script_path)
  }
}

#test2 <- create_counts()

# Merging all the data and creating the results file


merge_counts <- function() {
  # Protein coding genes
  rdata_protein_coding_genes <- load('data/protein_coding_genes.RData')
  protein_coding_genes <- get("protein_coding_genes") %>%
    dplyr::select(hgnc_id)

  # Input genes
  rdata_input_genes <- load('data/input_genes.RData')
  input_genes <- get("input_genes") %>%
    dplyr::select(hgnc_id)

  # Protein Complex
  rdata_protein_complexes <- load('data/complexportal_gene_counts.RData')
  protein_complexes <- get("complexportal_gene_counts") %>%
    dplyr::select(hgnc_id, ratio_input_genes_in_complexes)

  # Protein Families
  rdata_protein_families <- load('data/uniprot_pantherdb_gene_counts.RData')
  protein_families <- get("uniprot_pantherdb_gene_counts") %>%
    dplyr::select(hgnc_id, ratio_input_genes_in_families)

  # Pathways
  rdata_pathways <- load('data/reactome_input_genes_pathways_gene_counts.RData')
  pathways <- get("reactome_input_genes_pathways_gene_counts") %>%
    dplyr::select(hgnc_id, ratio_input_genes_in_pathways)

  # Paralogues
  rdata_paralogues <- load('data/biomart_paralogues_ratio.RData')
  paralogues <- get("biomart_paralogues_ratio") %>%
    dplyr::select(hgnc_id, ratio_paraloginputgenes_to_paralogs)

  # PPI
  rdata_ppi <- load('data/string_interactions_ratio.RData')
  ppi <- get("string_interactions_ratio") %>%
    dplyr::select(hgnc_id, ratio_interactioninputgenes_to_interactions)

  # Merging everything
  list_of_dfs <- list(protein_coding_genes, protein_complexes, protein_families,
                      pathways, paralogues, ppi)

  results <- reduce(list_of_dfs, left_join, by = 'hgnc_id') %>%
    dplyr::rename(protein_complex_score = ratio_input_genes_in_complexes) %>%
    dplyr::rename(protein_family_score = ratio_input_genes_in_families) %>%
    dplyr::rename(pathway_score = ratio_input_genes_in_pathways) %>%
    dplyr::rename(paralogue_score = ratio_paraloginputgenes_to_paralogs) %>%
    dplyr::rename(ppi_score = ratio_interactioninputgenes_to_interactions) %>%
    arrange(desc(protein_complex_score), desc(protein_family_score),
            desc(pathway_score), desc(paralogue_score), desc(ppi_score))
  #print(results)

  write.csv(results_with_previous_classes, file = 'results.csv')
}



# Comparing the results to our previous labels ---------------------------------

#previous_gene_classes <- data.frame(read_delim('./data/input/ndd_dataset_gene_classes.txt',
 #                                              delim = '\t')) %>%
#  dplyr::select(hgnc_id, ndd_ad_classes)
#print(previous_gene_classes)

#results_with_previous_classes <- results %>%
 # left_join(previous_gene_classes, by = 'hgnc_id')
#print(results_with_previous_classes)

# Removing duplicate rows
#results_with_previous_classes <- distinct(results_with_previous_classes)

# Saving the results -----------------------------------------------------------

#save_data_with_datetime('./data/results', 'protein_coding_genes_input_genes_5p_results',
 #                       results_with_previous_classes)





# Run the plots?







