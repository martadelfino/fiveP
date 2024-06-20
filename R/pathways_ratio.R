################################################################################
###### Script for creating counts of the pathways ##############################
################################################################################

library('tidyverse')

#source("./src/functions.R")

# Accessing the latest Reactome enriched pathways file -------------------------

#enriched_pathways <- access_last_saved_file('./data/raw/', "reactome_enriched")


# Accessing the latest Reactome pathways of input genes file -------------------

#input_genes_pathways <- access_last_saved_file('./data/processed/', "reactome_gene_list_Uniprot2Reactome")

rdata_input_genes_pathways <- load('data/reactome_gene_list_Uniprot2Reactome.RData')
input_genes_pathways <- get("reactome_gene_list_Uniprot2Reactome")


# Accessing the Uniprot2Reactome file of all genes in the pathways -------------

#Uniprot2Reactome <- access_last_saved_file('./data/processed/', "reactome_Uniprot2Reactome_cleaned")

rdata_Uniprot2Reactome <- load('data/reactome_Uniprot2Reactome_cleaned.RData')
Uniprot2Reactome <- get("reactome_Uniprot2Reactome_cleaned")


# Combining the enrichment and Uniprot2Reactome --------------------------------

#enriched_Uniprot2Reactome <- Uniprot2Reactome %>%
#  dplyr::filter(pathway_id %in% enriched_pathways$ID)
#print(enriched_Uniprot2Reactome)

# There are some proteins that don't map to HGNCs, and this is because they
# are immunology related. such as subsections of antibodies

# Saving the enriched pathways -------------------------------------------------

#save_data_with_datetime('./data/processed', 'reactome_enriched_pathways_',
#                        enriched_Uniprot2Reactome)


# Counting the number of input genes per pathway -------------------------------

#input_genes <- read_delim('./data/input/gene_list.txt', delim = '\t')
#print(input_genes)

rdata_input_genes <- load('data/input_genes.RData')
input_genes <- get("input_genes") %>%
  dplyr::select(hgnc_id)


#enriched_reactome_counts <- enriched_Uniprot2Reactome %>%
#  mutate(input_gene_yes_or_no = ifelse(hgnc_id %in% input_genes$hgnc_id, 1, 0))
#print(enriched_reactome_counts, n=5)

# Counting the number of input genes per pathway

#enriched_reactome_counts <- enriched_reactome_counts %>%
#  group_by(pathway_id) %>%
#  mutate(num_genes_in_pathway = n(),
#         numb_input_gene_per_pathway = sum(input_gene_yes_or_no))
#print(enriched_reactome_counts, n=5)

# Saving the enriched reactome pathways with number in input genes per pathway -

#save_data_with_datetime('./data/counts', 'reactome_enriched_pathways_inputgene_yes_or_no',
#                        enriched_reactome_counts)


# Creating a df for the number of unique genes in each pathway -----------------

# Counting the number of unique genes in each pathway that gene is related to
#enriched_reactome_counts_per_gene <- enriched_reactome_counts %>%
#  group_by(hgnc_id) %>%
#  dplyr::mutate(
#    num_pathways = n_distinct(pathway_id),
#    num_unique_genes_in_pathways = sum(length(unique(enriched_reactome_counts$hgnc_id[enriched_reactome_counts$pathway_id %in% pathway_id])) - 1),
#    num_input_genes_in_pathways = sum(unique(enriched_reactome_counts$hgnc_id[enriched_reactome_counts$pathway_id %in% pathway_id]) %in% input_genes$hgnc_id) - (hgnc_id %in% input_genes$hgnc_id)
#  ) %>%
#  dplyr::mutate(ratio_input_genes_in_pathways = num_input_genes_in_pathways / num_unique_genes_in_pathways) %>%
#  dplyr::select(hgnc_id, uniprot_ids, pathway_id, num_pathways,
#              num_unique_genes_in_pathways, num_input_genes_in_pathways,
#              ratio_input_genes_in_pathways) %>% arrange(hgnc_id)
#print(enriched_reactome_counts_per_gene, n= 10)

#enriched_reactome_counts_per_gene_final <- enriched_reactome_counts_per_gene %>%
#  dplyr::select(hgnc_id, uniprot_ids, num_pathways,
#                num_unique_genes_in_pathways, num_input_genes_in_pathways,
#                ratio_input_genes_in_pathways) %>% unique()
#print(enriched_reactome_counts_per_gene_final, n= 10)

# Saving the enriched reactome pathways with counts ----------------------------

#save_data_with_datetime('./data/counts', 'reactome_enrichedpathways_gene_counts',
#                        enriched_reactome_counts_per_gene_final)




# Now doing the same but with all the input gene pathways, not just enriched ones -----

# Combining the input genes pathways and Uniprot2Reactome ----------------------

input_genes_Uniprot2Reactome <- Uniprot2Reactome %>%
  dplyr::filter(pathway_id %in% input_genes_pathways$pathway_id)
print(input_genes_Uniprot2Reactome)

# Saving the input genes reactome pathways -------------------------------------

#save_data_with_datetime('./data/processed', 'reactome_input_genes_pathways_',
 #                       input_genes_Uniprot2Reactome)

save(input_genes_Uniprot2Reactome, file = "data/reactome_input_genes_pathways.RData")


# Counting the number of input genes per pathway -------------------------------

# Checking if the individual genes are input genes or not
reactome_counts <- input_genes_Uniprot2Reactome %>%
  mutate(input_gene_yes_or_no = ifelse(hgnc_id %in% input_genes$hgnc_id, 1, 0))
print(reactome_counts, n=5)

# Counting the number of input genes per pathway
reactome_counts <- reactome_counts %>%
  group_by(pathway_id) %>%
  mutate(num_genes_in_pathway = n(),
         numb_input_gene_per_pathway = sum(input_gene_yes_or_no))
print(reactome_counts, n=5)


# Saving the enriched reactome pathways with number in input genes per pathway  ------------------------------------------------------------------------

#save_data_with_datetime('./data/counts/', 'reactome_input_genes_inputgene_yes_or_no',
 #                       reactome_counts)

save(reactome_counts, file = "data/reactome_input_genes_inputgene_yes_or_no.RData")


# Creating a df for the number of unique genes in each pathway -----------------

# Counting the number of unique genes in each pathway that gene is related to
reactome_counts_per_gene <- reactome_counts %>%
  group_by(hgnc_id) %>%
  dplyr::mutate(
    num_pathways = n_distinct(pathway_id),
    num_unique_genes_in_pathways = sum(length(unique(reactome_counts$hgnc_id[reactome_counts$pathway_id %in% pathway_id])) - 1),
    num_input_genes_in_pathways = sum(unique(reactome_counts$hgnc_id[reactome_counts$pathway_id %in% pathway_id]) %in% input_genes$hgnc_id) - (hgnc_id %in% input_genes$hgnc_id)
  ) %>%
  dplyr::mutate(ratio_input_genes_in_pathways = num_input_genes_in_pathways / num_unique_genes_in_pathways) %>%
  dplyr::select(hgnc_id, uniprot_ids, pathway_id, num_pathways,
                num_unique_genes_in_pathways, num_input_genes_in_pathways,
                ratio_input_genes_in_pathways) %>% arrange(hgnc_id)
print(reactome_counts_per_gene, n= 10)

reactome_counts_per_gene_final <- reactome_counts_per_gene %>%
  dplyr::select(hgnc_id, uniprot_ids, num_pathways,
                num_unique_genes_in_pathways, num_input_genes_in_pathways,
                ratio_input_genes_in_pathways) %>% unique() %>%
  dplyr::mutate(ratio_input_genes_in_pathways = ifelse(is.na(ratio_input_genes_in_pathways),
                                                       0, ratio_input_genes_in_pathways))
print(reactome_counts_per_gene_final, n= 1000)

#duplicate_rows <- data.frame(reactome_counts_per_gene_final) %>%
#  filter(duplicated(.))
#print(duplicate_rows)

# Saving the pathway counts per gene -------------------------------------------

#save_data_with_datetime('./data/counts', 'reactome_input_genes_pathways_gene_counts',
 #                       reactome_counts_per_gene_final)

save(reactome_counts_per_gene_final,
     file = "data/reactome_input_genes_pathways_gene_counts.RData")




# Plots ------------------------------------------------------------------------


# make plots to compare the results of the enriched pathways and all the pathways

