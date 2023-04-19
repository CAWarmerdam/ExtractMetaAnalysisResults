#!/usr/bin/env Rscript

## ----
## Author:  C.A. (Robert) Warmerdam
## Email:   c.a.warmerdam@umcg.nl
##
## Copyright (c) C.A. Warmerdam, 2021
##
## This program is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.
##
## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.
##
## A copy of the GNU General Public License can be found in the LICENSE file in the
## root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
## ----

# Load libraries
library(argparse)
library(tidyverse)
library(DBI)
library(arrow)
#library(duckdb)

# Declare constants

parser <- ArgumentParser(description='')
parser$add_argument('--dataset',
                    help='Path to dataset')
parser$add_argument('--variants',
                    help='Variants to sample')

gene_selection <- c(
  "ENSG00000090339")

# Declare function definitions

# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  # Process input
  args <- parser$parse_args(argv)
  eqtls_dataset <- arrow::open_dataset(args$dataset)

  to.int <- function(values, data.type) {
    return(unlist(sapply(values, function(value) Scalar$create(value, data.type))))
  }

  variant_tibble <- read_delim(args$variants, delim=" ") %>% 
    transmute(chromosome = CHR, bp = bp, variant=ID) %>%
    as_arrow_table(
      schema = schema(chromosome = int8(), bp = int32(), variant=string()))

  print(variant_tibble)
  
  # Select results for specific phenotype and time this.
  for (gene in gene_selection) {
  
    start_time <- Sys.time()
    selection <- eqtls_dataset %>% filter(phenotype == gene) %>% collect() %>%
      mutate(p = pnorm(abs(beta/standard_error), lower.tail=FALSE)*2)
    end_time <- Sys.time()
    arrow_time <- end_time - start_time
    message(sprintf("Filtering gene: %s", format(arrow_time)))

    write.table(selection, sprintf("ExtractedPhenotype_%s.txt.gz", gene), row.names=F, sep="\t", quote=F)
 
  }

  #start_time <- Sys.time()
  #selection_positional <- eqtls_dataset %>%
  #  semi_join(variant_tibble, by = c("chromosome"="chromosome", "bp"="bp")) %>% collect()
  #end_time <- Sys.time()
  #arrow_time <- end_time - start_time
  #message(sprintf("Filtering gene: %s", format(arrow_time)))

  #
  # # Perform method
  # selection <- DBI::dbGetQuery(
  #   eqtls_db_connection,
  #   "SELECT * FROM eqtls")

  # Process output
  #write.table(selection, "output.txt")
  #write.table(selection_positional, "output_positional.txt")
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}

