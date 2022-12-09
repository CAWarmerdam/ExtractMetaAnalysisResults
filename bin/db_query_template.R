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

  # Select results for specific phenotype and time this.
  start_time <- Sys.time()
  selection <- eqtls_dataset %>% filter(phenotype == "ENSG00000164167") %>% head() %>% collect()
  end_time <- Sys.time()
  arrow_time <- end_time - start_time
  message(sprintf("Filtering gene: %s", format(arrow_time)))

  #
  # # Perform method
  # selection <- DBI::dbGetQuery(
  #   eqtls_db_connection,
  #   "SELECT * FROM eqtls")

  # Process output
  write.table(selection, "output.txt")
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
