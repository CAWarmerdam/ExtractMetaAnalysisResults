#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(tidyverse)
library(dtplyr)
library(arrow)

# Declare constants


# Declare function definitions
extract_locus_as_dt <- function(permuted_dataset_path, variant_index_end, variant_index_start) {
  ds <- arrow::open_dataset(permuted_dataset_path)
  z_score_dt <- ds %>%
    filter(between(variant_index, variant_index_start, variant_index_end)) %>%
    mutate(z_score = beta / standard_error) %>%
    as.data.table()
  return(z_score_dt)
}


# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input
  variant_ref <- "/1000G-30x_index.parquet"
  position_chromosome <- 21
  position_start <- 34446503
  position_end <- 37727357

  # Get the variants to load
  locus_variant_reference <- arrow::read_parquet(variant_ref) %>%
    filter(chromosome == position_chromosome, between(bp, position_start, position_end)) %>% collect()

  # Get the indices of the variants to laod
  variant_index_start <- min(locus_variant_reference$variant_index)
  variant_index_end <- max(locus_variant_reference$variant_index)

  # Get extract data for locus to calculate LD for.
  start.time.shared <- Sys.time() # Record start time

  # Open permuted parquet dataset from shared filesystem
  permuted_dataset_path <- "/meta"
  z_score_dt_shared <- extract_locus_as_dt(permuted_dataset_path, variant_index_end, variant_index_start)

  end.time.shared <- Sys.time() # Record end time

  time.taken.shared <- end.time.shared - start.time.shared # Calculate time taken
  print("Duration using shared filesystem:")
  print(time.taken.shared)

  # Get extract data for locus to calculate LD for.
  start.time.local <- Sys.time() # Record start time

  # Open locally synced permuted parquet dataset
  permuted_dataset_path <- "/local_meta"
  z_score_dt_local <- extract_locus_as_dt(permuted_dataset_path, variant_index_end, variant_index_start)

  end.time.local <- Sys.time() # Record end time

  time.taken.local <- end.time.local - start.time.local # Calculate time taken
  print("Duration for local copy:")
  print(time.taken.local)

  # Get z-score matrix from data table
  z_score_mat <- z_score_dt_local %>%
    pivot_wider(id_cols = "variant_index", names_from = "phenotype", values_from = "z_score") %>%
    collect() %>% as.data.table() %>% as.matrix(rownames=1)

  start.time <- Sys.time()

  z_score_mat <- z_score_mat - rowMeans(z_score_mat)
  # Standardize each variable
  z_score_mat <- z_score_mat / sqrt(rowSums(z_score_mat^2))
  # Calculate correlations
  ld_matrix <- tcrossprod(z_score_mat)

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  print(time.taken)


  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}