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
library(data.table)
library(tidyverse)

# Declare constants

# Declare function definitions

# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  # Load uncorrelated variants, MAF data
  frequencies <- as_tibble(fread(argv[1], data.table=F))
  vars_ref <- as_tibble(fread(argv[2], data.table=F, header=F))
  results_target <- as_tibble(fread(argv[3], data.table=F, header=F))

  print(head(frequencies))

  # Calculate MAF
  freq_tib_maf <- frequencies %>% mutate(MAF = pmin(AF, 1-AF))

  print(head(freq_tib_maf))

  print(min(freq_tib_maf$MAF))

  # Filter frequencies according to the independent SNPs.
  quantile_results <- quantile(freq_tib_maf %>% filter(freq_tib_maf$ID %in% vars_ref$V1) %>% pull(MAF), 1:10/10)

  print(quantile_results)
  write.table(unname(quantile_results), paste0(argv[4], ".breakpoints.txt"), quote=F, col.names =F, row.names=F)

  freq_tib_binned <- freq_tib_maf %>%
    filter(ID %in% results_target$V1) %>%
    mutate(bin = cut(MAF, breaks=c(0.01, quantile_results), labels=FALSE),
           lower=c(0.01, quantile_results)[bin],
           upper=quantile_results[bin])

  for (i in seq_along(quantile_results)) {
    write.table(freq_tib_binned %>% filter(bin==i) %>% pull(ID),
                paste0(argv[4], ".var_set_bin_", i, ".txt"), col.names=F, row.names=F, quote=F)
  }
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}

