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
library(susieR)

# Declare constants

# Create a parser object
parser <- ArgumentParser(description = 'Run SuSiE over loci.')

# Add command-line arguments
parser$add_argument(
  "--permuted",
  required = TRUE,
  help = "Path to the directory containing permuted data."
)

parser$add_argument(
  "--empirical",
  required = TRUE,
  help = "Path to the directory containing empirical data."
)

parser$add_argument(
  "--variant-reference",
  required = TRUE,
  help = "Path to the variant reference file."
)

parser$add_argument(
  "--uncorrelated-genes",
  required = TRUE,
  help = "File containing uncorrelated genes."
)

parser$add_argument(
  "--bed-files",
  required = TRUE,
  nargs = "+",  # Allows multiple BED files to be passed as arguments
  help = "Space-separated list of BED files."
)

# Declare function definitions

#' Convert Z-Score to Correlation Coefficient
#'
#' This function converts a Z-score back to a correlation coefficient `r` based on
#' the given degrees of freedom `df`. This is useful for comparing Z-scores when sample sizes
#' are dissimilar
#'
#' @param z A numeric vector of Z-scores.
#' @param df A numeric vector of degrees of freedom corresponding to the Z-scores.
#' @return A numeric vector of correlation coefficients. If the length of `z` and `df` do not match, the function returns 0.
#' @details
#' The function calculates the t-value from the Z-score using the inverse cumulative distribution
#' function of the normal distribution. It then converts the t-value to the corresponding
#' correlation coefficient. The sign of the correlation coefficient is adjusted based on the
#' sign of the original Z-score.
#'
#' @examples
#' # Convert a Z-score of 2.5 with 10 degrees of freedom to a correlation coefficient
#' zToCor(2.5, 10)
#'
#' # Convert multiple Z-scores to correlation coefficients
#' zToCor(c(1.5, -2.0), c(15, 20))
#'
#' @export
zToCor <- function(z, df){
  if(length(z) != length(df)){
    return(0L)  #Error condition
  }
  t <- qt(pnorm(-abs(z),log.p = T), df, log.p=T)
  r <- t/sqrt(df+t^2)
  r <- ifelse(z > 0, r * -1, r)
  return(r)
}

# Extract locus as data table
get_ld_matrix <- function(permuted_dataset, variant_index_start, variant_index_end) {
  z_score_dt <- permuted_dataset %>%
    filter(between(variant_index, variant_index_start, variant_index_end)) %>%
    mutate(z_score = beta / standard_error) %>%
    as.data.table()

  print(as_tibble(z_score_dt))

  # Get z-score matrix from data table
  z_score_mat <- z_score_dt %>%
    pivot_wider(id_cols = "variant_index", names_from = "phenotype", values_from = "z_score") %>%
    collect() %>% as.data.table() %>% as.matrix(rownames=1)

  # Get z-score matrix from data table
  n_mat <- z_score_dt %>%
    pivot_wider(id_cols = "variant_index", names_from = "phenotype", values_from = "sample_size") %>%
    collect() %>% as.data.table() %>% as.matrix(rownames=1)

  # Convert Z-scores to pearson correlations, using sample size - 1 as the degrees of freedom
  rho_mat <- zToCor(z_score_mat, n_mat-1)

  start.time <- Sys.time()

  rho_mat <- rho_mat - rowMeans(rho_mat)
  # Standardize each variable
  rho_mat <- rho_mat / sqrt(rowSums(rho_mat^2))
  # Calculate correlations
  ld_matrix <- tcrossprod(rho_mat)

  end.time <- Sys.time()
  time.taken <- end.time - start.time

  print(time.taken)
  return(ld_matrix)
}

finemap_locus <- function(empirical_dataset, permuted_dataset, locus_bed, variant_reference) {
  locus_chromosome <- unique(locus_bed %>% pull(chromosome))
  locus_start <- min(locus_bed %>% pull(start))
  locus_end <- max(locus_bed %>% pull(end))

  # Get the variants to load
  locus_variant_reference <- variant_reference %>%
    filter(chromosome == locus_chromosome, between(bp, locus_start, locus_end)) %>% collect()

  # Get the indices of the variants to laod
  variant_index_start <- min(locus_variant_reference$variant_index)
  variant_index_end <- max(locus_variant_reference$variant_index)

  message("Starting to calculate LD...")
  start.time <- Sys.time()
  ld_matrix <- get_ld_matrix(permuted_dataset, variant_index_start, variant_index_end)
  print(as_tibble(ld_matrix))
  print(ld_matrix[1:10, 1:10])
  variant_order <- rownames(ld_matrix)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message("LD calculation done!")
  print(time.taken)

  # Do the remaining bit for finemapping the locus
  fine_mapping_results <- mapply(function(gene, start, end) {

    # Get the variants to load
    gene_locus_variant_reference <- variant_reference %>%
      filter(chromosome == locus_chromosome, between(bp, start, end)) %>% collect()

    # Get the indices of the variants to laod
    gene_variant_index_start <- min(gene_locus_variant_reference$variant_index)
    gene_variant_index_end <- max(gene_locus_variant_reference$variant_index)

    gene_summary_stats <- empirical_dataset %>%
      filter(phenotype == gene, between(variant_index, gene_variant_index_start, gene_variant_index_end)) %>%
      as.data.table()
    print(as_tibble(gene_summary_stats))

    # Yet to order the gene_summary_stats so that the variant ordering matches that of the ld_matrix
    # Do fine-mapping
    return(data.frame(list("gene" = c(), "variant" = c(), "cs" = c(), "pip" = c())))
  }, locus_bed$gene, locus_bed$start, locus_bed$end, SIMPLIFY =F)

  # Return
  return(fine_mapping_results)
}


# Main
print = T
direct = T
formatted = T

if(formatted){
  outMat = NULL
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

  # Access the arguments
  cat("Permuted data directory:", args$permuted, "\n")
  cat("Empirical data directory:", args$empirical, "\n")
  cat("Variant reference file:", args$variant_reference, "\n")
  cat("Uncorrelated genes file:", args$uncorrelated_genes, "\n")
  cat("BED files:", paste(args$bed_files, collapse = ", "), "\n")

  # get path of parquet from arguments
  sumstats_path <- args$empirical
  permuted_dataset_path <- args$permuted
  variant_reference_path <- args$variant_reference

  variant_reference <- arrow::read_parquet(variant_reference_path)

  # load in parquet with Robert's strategy
  empirical_dataset <- open_dataset(sumstats_path)
  permuted_dataset <- arrow::open_dataset(permuted_dataset_path)

  fine_mapping_results_per_locus <- mapply(function(bed_file) {
    locus_bed <- fread(bed_file, col.names = c("chromosome", "start", "end", "gene", "cluster"))
    # Assumes fine_mapping_output is some sort of data.frame/data.table or tibble
    fine_mapping_output <- finemap_locus(empirical_dataset=empirical_dataset,
                                         permuted_dataset=permuted_dataset,
                                         locus_bed=locus_bed,
                                         variant_reference=variant_reference)
    return(fine_mapping_output)
  }, args$bed_files, SIMPLIFY = F)

  combined_results <- bind_rows(fine_mapping_results_per_locus)

  fwrite(combined_results, "fine_mapping_results.tsv", sep="\t", quote=F, row.names=F, col.names=T)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
