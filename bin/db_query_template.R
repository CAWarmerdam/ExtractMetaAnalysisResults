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

# Declare constants

parser <- ArgumentParser(description='')
parser$add_argument('--database',
                    help='Path to database')

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
  eqtls_db_connection <- DBI::dbConnect(RSQLite::SQLite(), dbname=args$database)

  DBI::dbExecute(eqtls_db_connection, "SELECT load_extension('/tools/libparquet.so');")

  # Perform method
  selection <- DBI::dbGetQuery(
    eqtls_db_connection,
    "SELECT * FROM eqtls")

  # Process output
  write.table(selection, "output.txt")
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
