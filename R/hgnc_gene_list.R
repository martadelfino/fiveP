################################################################################
###### Script to access protein coding genes list HGNC #########################
################################################################################

#if (!require("dplyr")) install.packages("dplyr")
library("tidyverse")

#if (!require("logr")) install.packages("logr")
library("logr")

# hgnc gene file ----------------------------------------------------------

# If the file doesn't already exist, we will read it via FTP
filename <- file.path("gene_with_protein_product.txt")

if (!file.exists(filename)) {
  filename <- paste("ftp://ftp.ebi.ac.uk/pub/databases/genenames/hgnc/tsv/locus_types",
                    "gene_with_protein_product.txt",
                    sep = "/")
}

# filepath <- system.file("data", "gene_with_protein_product.txt", package = "your_package_name")

protein_coding_genes <- readr::read_delim(filename,
                                          delim = "\t",
                                          col_names = TRUE) %>%
  as.data.frame()


# Save file ------------------------------------------------------------------

save(protein_coding_genes, file = "data/protein_coding_genes.RData")

#write.table(protein_coding_genes, "./data/protein_coding_genes.txt",
 #           quote = F, sep = "\t", row.names = F)

#write.table(protein_coding_genes, "./data/raw/protein_coding_genes.txt",
 #           quote = F, sep = "\t", row.names = F)
