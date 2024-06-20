################################################################################
###### Script for obtaining protein complex participants #######################
################################################################################

library('tidyverse')

#source("./src/functions.R")


# Accessing the latest Uniprot results file ------------------------------------

#uniprot_input_gene_symbol_results_cleaned <- access_last_saved_file('./data/processed/', 'uniprot_input_gene_symbol_results_cleaned')
#print(uniprot_input_gene_symbol_results_cleaned$ComplexPortal)

rdata_uniprot_input_gene_symbol_results_cleaned <- load('data/uniprot_input_gene_symbol_results_cleaned.RData')
uniprot_input_gene_symbol_results_cleaned <- get("uniprot_input_gene_symbol_results_cleaned")

check_for_weird_characters(uniprot_input_gene_symbol_results_cleaned$ComplexPortal)

# Creating a new df of the complexes -------------------------------------------

input_genes_protein_complexes_expanded <- uniprot_input_gene_symbol_results_cleaned %>%
  tidyr::separate_rows(ComplexPortal, sep = ";") %>% distinct() %>%
  filter(ComplexPortal != "") %>%
  dplyr::select(ComplexPortal) %>% distinct() %>%
  dplyr::rename(complex_id = ComplexPortal)

# removing extra bits
input_genes_protein_complexes_expanded$complex_id <- trimws(sub("\\[.*?\\]", "", input_genes_protein_complexes_expanded$complex_id))
print(input_genes_protein_complexes_expanded, n=100)


# Querying ComplexPortal -------------------------------------------------------

ComplexPortal <- read_delim('https://ftp.ebi.ac.uk/pub/databases/intact/complex/current/complextab/9606.tsv')
print(ComplexPortal$`#Complex ac`)
print(head(ComplexPortal$`Identifiers (and stoichiometry) of molecules in complex`))
print(head(ComplexPortal$`Expanded participant list`))
print(head(ComplexPortal$Comment))
print(head(ComplexPortal))

#test2 <- ComplexPortal[grepl('CPX-3195', ComplexPortal$`#Complex ac`), ]
#print(test2$`Identifiers (and stoichiometry) of molecules in complex`)

# Save the raw file ------------------------------------------------------------

#save_data_with_datetime('./data/raw/', 'complexportal', ComplexPortal)

save(ComplexPortal, file = "data/complexportal.RData")


# Cleaning ComplexPortal file --------------------------------------------------

# Note: the ComplexPortal is updated every month, so the date will be important

# Selecting and renaming required columns
ComplexPortal_participants <- ComplexPortal %>%
  dplyr::select(`#Complex ac`, `Expanded participant list`) %>%
  dplyr::rename(complex_id = `#Complex ac`) %>%
  dplyr::rename(uniprot_ids = `Expanded participant list`)
print(head(ComplexPortal_participants))

# Creating a new row for each protein
ComplexPortal_participants_separated <- ComplexPortal_participants %>%
  separate_rows(uniprot_ids, sep = '\\|')
print(ComplexPortal_participants_separated)

# Removing extra information enclosed in '[]' and after '-'
ComplexPortal_participants_separated$uniprot_ids <- gsub("\\[.*?\\]", "", ComplexPortal_participants_separated$uniprot_ids)
ComplexPortal_participants_separated$uniprot_ids <- sub("-.*", "", ComplexPortal_participants_separated$uniprot_ids)

# Removing the stochiometry information in the brackets
ComplexPortal_participants_separated$uniprot_ids <- trimws(sub("\\(.*?\\)", "", ComplexPortal_participants_separated$uniprot_ids))
print(ComplexPortal_participants_separated)

# Removing rows without a uniprot_id (these would've been other things like molecules)
ComplexPortal_participants_separated <- subset(ComplexPortal_participants_separated, uniprot_ids != "")

#test3 <- ComplexPortal_participants_separated[grepl('CPX-6904', ComplexPortal_participants_separated$complex_id), ]
#print(test3)

# Checking for any other weird characters
check_for_weird_characters(ComplexPortal_participants_separated$uniprot_ids)


# Only keeping the complexes identified in the input gene lists ----------------

# Joining the input genes protein complexes with the participants for each complex
input_genes_complexportal_participants <- input_genes_protein_complexes_expanded %>%
  left_join(ComplexPortal_participants_separated, by = 'complex_id' ,
  relationship = "many-to-many")

print(input_genes_complexportal_participants)


### Adding hgnc_id

# getting the protein coding genes file
#hgnc_uniprot_symbol <- read_delim('./data/raw/protein_coding_genes.txt', '\t') %>%
 # dplyr::select(hgnc_id, uniprot_ids, symbol)
rdata_hgnc_uniprot_symbol <- load('data/protein_coding_genes.RData')
hgnc_uniprot_symbol <- get("protein_coding_genes") %>%
  dplyr::select(hgnc_id, uniprot_ids, symbol)

hgnc_uniprot_symbol$uniprot_ids <- trimws(sub('\\|.*', "", hgnc_uniprot_symbol$uniprot_ids))

input_genes_complexportal_participants_hgnc <- input_genes_complexportal_participants %>%
  left_join(hgnc_uniprot_symbol, by = 'uniprot_ids',
            relationship = "many-to-many") %>% unique()
print(input_genes_complexportal_participants_hgnc)

rows_with_na <- input_genes_complexportal_participants_hgnc[!complete.cases(input_genes_complexportal_participants_hgnc), ]
print(rows_with_na)

# there are still 5 proteins with no HGNC matching. I could use the checker Pilar made?


# Saving the cleaning ComplexPortal file ---------------------------------------

#save_data_with_datetime('./data/processed/', 'complexportal_cleaned',
 #                       input_genes_complexportal_participants_hgnc)

save(input_genes_complexportal_participants_hgnc,
     file = "data/complexportal_cleaned.RData")


# Calculating how many complexes are input gene list complexes from the total complexes ----

total_complexes <- ComplexPortal %>%
  dplyr::select('#Complex ac') %>% rename(complex_id = '#Complex ac') %>%
  distinct()
# there are 2050 total protein complexes

complexes_that_are_input_gene_complexes <- input_genes_complexportal_participants %>%
  dplyr::select('complex_id') %>% distinct()
# 343 are NDD complexes



