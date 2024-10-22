# Getting the initial 5P data objects and downloading them as Rdata for Pilar


library(fiveP)
library(tidyverse)

gene_classes <- readr::read_delim('CodeReview_data.txt', delim = '\t',
                                  show_col_types = FALSE)

input_genes <- gene_classes %>%
  dplyr::filter(ndd_ar_classes == 'positive') %>%
  dplyr::select(hgnc_id)

hgnc_gene_list <- fetch_hgnc_gene_list()

paralogues <- fetch_paralogues(hgnc_gene_list)
pathways <- fetch_pathways(hgnc_gene_list, input_genes)
ppi <- fetch_ppi(hgnc_gene_list)
uniprot <- fetch_uniprot(hgnc_gene_list, input_genes)
protein_complex <- fetch_protein_complex(hgnc_gene_list, uniprot)
protein_families <- fetch_protein_families(hgnc_gene_list, uniprot)


save(input_genes, file = "input_genes.RData")
save(hgnc_gene_list, file = "hgnc_gene_list.RData")
save(paralogues, file = "paralogues.RData")
save(pathways, file = "pathways.RData")
save(ppi, file = "ppi.RData")
save(uniprot, file = "uniprot.RData")
save(protein_complex, file = "protein_complex.RData")
save(protein_families, file = "protein_families.RData")
