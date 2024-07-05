# fiveP

## Overview

The main function, `get_fiveP(input_genes)`, takes a df of HGNC IDs. The function outputs a dataframe with HGNC IDs for all protein coding genes, with five columns (paralogues, pathways, PPI, protein complexes, protein families) of scores. The scores are calculated based on how similar each protein coding gene is to the input genes provided by the user.

## Installation 

```{r}
install.packages("devtools")

devtools::install_github("martadelfino/fiveP")
```

## Usage 

Example

```{r}
library(fiveP)
library(tidyverse)

input_genes <- tibble(hgnc_id = c("HGNC:19743", 
                                  "HGNC:9202",
                                  "HGNC:8653",
                                  "HGNC:6936",
                                  "HGNC:4878"))

result <- get_fiveP(input_genes)
head(result)
```

-   More examples in `biallelic.Rmd`

## Features 

-   Main function: `get_fiveP()`

## Documentation 

-   Roxygen docs in `man/`
