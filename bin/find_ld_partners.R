#!/usr/bin/env Rscript

## ----
## Author:  T. van Lieshout
## Email:   t.van.lieshout@umcg.nl
##
## Copyright (c) T. van Lieshout, 2024
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
library(data.table)
library(tidyverse)
library(dtplyr)
library(arrow)

# Declare constants

# Create a parser object
parser <- ArgumentParser(description = 'Run SuSiE over loci.')

parser$add_argument(
  "--independent-variants",
  required=T,
  help = "Path to putatively independent variants"
)

parser$add_argument(
  "--variant-reference",
  required = TRUE,
  help = "Path to the variant reference file."
)

# Add command-line arguments
parser$add_argument(
  "--ld",
  required = TRUE,
  help = "Path to the directory containing ld data."
)

parser$add_argument(
  "--ld-type",
  required = FALSE,
  default="pcs",
  choices=c("gene-set", "pcs", "dosages", "feather"),
  help = "Type of LD"
)

parser$add_argument(
  "--r2-threshold",
  required = TRUE,
  default=0.9,
  type='numeric'
)


# Declare function definitions

ld_calculator <- function(rho_mat) {
  rho_mat <- rho_mat - rowMeans(rho_mat)
  # Standardize each variable
  rho_mat <- rho_mat / sqrt(rowSums(rho_mat^2))
  # Calculate correlations
  ld_matrix <- tcrossprod(rho_mat)
  return(ld_matrix)
}


susie_partners <- function(lead_variant_index, lead_variant_chromosome, lead_variant_bp, variant_reference, ld_func, r2_threshold, kb_threshold) {
  variants_in_window <- variant_reference %>%
    filter(chromosome == lead_variant_chromosome,
           between(bp, lead_variant_bp-kb_threshold, lead_variant_bp+kb_threshold)) %>% collect()

  cat(lead_variant_chromosome, lead_variant_bp)

  ld_matrix <- ld_func(
    variants_in_window$chromosome[1],
    min(variants_in_window$variant_index),
    max(variants_in_window$variant_index))

  ld_matrix <- ld_matrix[rownames(ld_matrix), rownames(ld_matrix)]

  ld_values <- tibble(
    variant_index_tagging = as.integer(rownames(ld_matrix)),
    r_tagging = ld_matrix[as.character(lead_variant_index), ]) %>%
    filter(r_tagging^2 > r2_threshold)

  return(list(
    "variant_index_tagging" = ld_values$variant_index_tagging,
    "r_tagging" = ld_values$r_tagging))
}


#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  # Parse the arguments
  args <- parser$parse_args()

  ld_path <- "/tmp/ld_panel_10497_pcs/"
  ld_type <- "pcs"
  susie_out_path <- "independent_variants_filtered_lbf2_mlog10p5_annotated_20250325.txt.gz"
  min_r2_threshold <- 0.9
  variant_reference_path <- "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet"

  # Access the arguments
  ld_path <- args$ld
  cat("LD data directory:", ld_path, "\n")
  ld_type <- args$ld_type
  cat("LD type:", ld_type, "\n")
  susie_out_path <- args$independent_variants
  cat("Input table:", susie_out_path, "\n")
  min_r2_threshold <- args$r2_threshold
  cat("R2 threshold", min_r2_threshold, "\n")
  cat("Variant reference":variant_reference_path)

  variant_reference_path <- args$variant_reference

  # get path of parquet from arguments
  variant_reference <- arrow::read_parquet(variant_reference_path)
  variant_table <- fread(susie_out_path) %>%
    inner_join(variant_reference)

  permuted_dataset_path <- ld_path
  min_r2_threshold

  # Switched to new ld reference dataset files with just phenotypes, variant_indices, and rho values
  # Added filtering on uncorrelated genes, since the dataset is no longer prefiltered on uncorrelated genes
  # (the entire set of permuted genes is in this dataset)
  permuted_dataset <- arrow::open_dataset(permuted_dataset_path)
  ld_func_window <- function(chromosome, variant_index_start, variant_index_stop) {
    rho_mat <- permuted_dataset %>% filter(
      chr == chromosome,
      between(variant_index, variant_index_start, variant_index_stop)) %>%
      collect() %>% as.data.table() %>%
      as.matrix(rownames=1)

    return(ld_calculator(rho_mat))
  }

  ld_partner_table <- variant_table %>%
    rowwise() %>%
    mutate(ld_partners = list(susie_partners(
      lead_variant_index = variant_index, lead_variant_chromosome = chromosome, lead_variant_bp = bp, variant_reference = variant_reference,
      ld_func = ld_func_window, r2_threshold = min_r2_threshold, kb_threshold = 0.25e6))) %>%
    unnest_wider(ld_partners) %>% unnest_longer(c(variant_index_tagging, r_tagging))

  fwrite(ld_partner_table, sprintf("%s_ld_partners.txt.gz", susie_out_path), sep="\t", col.names=T, row.names=F, quote=F)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
