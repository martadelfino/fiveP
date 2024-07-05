#if (!require("dplyr")) install.packages("dplyr")
library("tidyverse")

#if (!require("logr")) install.packages("logr")
library("logr")

#' Fetch all protein coding genes data
#'
#' This function fetches data from EBI to get all protein coding genes.
#'
#' @return A dataframe with protein coding gene data from EBI database.
#' @export
fetch_hgnc_gene_list <- function() {

  # Fetch hgnc gene file
  protein_coding_genes <- readr::read_delim("ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types/gene_with_protein_product.txt",
                                            delim = "\t",
                                            col_names = TRUE) %>%
    as.data.frame()

  cat('\n(1/12) finished running hgnc_gene_list.R\n')
  return(protein_coding_genes)

}
