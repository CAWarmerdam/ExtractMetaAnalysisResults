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
  "--ld",
  required = TRUE,
  help = "Path to the directory containing ld data."
)

parser$add_argument(
  "--ld-type",
  required = TRUE,
  choices=c("gene-set", "pcs", "dosages", "feather"),
  help = "Type of LD"
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
  required = FALSE,
  help = "File containing uncorrelated genes."
)

parser$add_argument(
  "--max-i2",
  required = FALSE, default=40, type = 'double',
  help = "Maximum i2"
)

parser$add_argument(
  "--min-n-prop",
  required = FALSE, default=0.8, type = 'double',
  help = "Minimum sample size proportion compared to maximum in locus"
)

parser$add_argument(
  "--no-adjust-stats",
  required = FALSE,
  default = FALSE,
  action = 'store_true',
  help = "Flag to disable adjusting betas and standard errors for sample size"
)

parser$add_argument(
  "--bed-files",
  required = TRUE,
  nargs = "+",  # Allows multiple BED files to be passed as arguments
  help = "Space-separated list of BED files."
)

parser$add_argument(
  "--debug",
  required = FALSE,
  default = FALSE,
  help = "Enables writing SuSIE input, choose from 'RSparsePro', 'locus-wise'"
)

parser$add_argument(
  "--dry-run",
  required = FALSE,
  default = FALSE,
  action = 'store_true',
  help = "Runs code without actually running SuSIE"
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
get_ld_matrix_alt <- function(permuted_dataset, variant_index_start, variant_index_end) {
  rho_mat <- permuted_dataset %>%
    filter(between(variant_index, variant_index_start, variant_index_end)) %>%
    collect() %>% as.data.table() %>% dcast(variant_index ~ phenotype, value.var = "rho") %>% 
    as.matrix(rownames=1)

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

# Extract locus as data table
get_ld_matrix_wide <- function(permuted_dataset, variant_index_start, variant_index_end) {
  rho_mat <- permuted_dataset %>%
    filter(between(variant_index, variant_index_start, variant_index_end)) %>%
    collect() %>% as.data.table() %>%
    as.matrix(rownames=1)

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

extract_summary_statistics <- function(gene_cluster, empirical_dataset, variant_reference) {

  cluster_summary_stats_list <- list()
  gene <- gene_cluster$gene[1]
  print(gene)

  message(sprintf("Extracting summary statistics for gene %s, %s loci", gene, nrow(gene_cluster)))
  for (i in seq_len(nrow(gene_cluster))) {
    locus_chromosome <- gene_cluster$chromosome[i]
    start <- gene_cluster$start[i]
    end <- gene_cluster$end[i]
    message(sprintf("%s:%s-%s", locus_chromosome, start, end))

    # Get variants in gene region
    gene_locus_variant_reference <- variant_reference %>%
      filter(chromosome == locus_chromosome,
             bp >= start,
             bp <= end)

    gene_variant_index_start <- min(gene_locus_variant_reference$variant_index, na.rm = TRUE)
    gene_variant_index_end <- max(gene_locus_variant_reference$variant_index, na.rm = TRUE)
    print(gene_variant_index_start)
    print(gene_variant_index_end)
    print(str(gene))
    print(empirical_dataset)

    # Read Parquet data with filtering
    query_result <- empirical_dataset %>%
      filter(phenotype == gene,
             between(variant_index, gene_variant_index_start, gene_variant_index_end)) %>%
      collect() %>% as.data.table()
    print(query_result)
    cluster_summary_stats_list[[as.character(i)]] <- query_result
    message(sprintf("Found %s variants", nrow(cluster_summary_stats_list[[as.character(i)]])))
  }
  summary_stats_bound <- bind_rows(cluster_summary_stats_list)
  message(sprintf("Total number of variants for %s: %s", gene, nrow(summary_stats_bound)))
  return(summary_stats_bound)
}

finemap_locus <- function(empirical_dataset, ld_func, locus_bed, variant_reference, min_sample_size_prop=0.8, max_i_squared=40, normalize_sumstats=T, debug=FALSE, dry_run=FALSE, nCS = 10) {
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
  ld_matrix <- ld_func(variant_index_start, variant_index_end)
  print(str(ld_matrix))
  variant_order <- rownames(ld_matrix)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message("LD calculation done!")
  print(time.taken)

  if (debug == 'locus-wise') {
    saveRDS(as.matrix(ld_matrix), sprintf("LD_%s:%d-%d.rds", locus_chromosome, locus_start, locus_end))
  }

  # Do the remaining bit for finemapping the locus
  fine_mapping_results <- bind_rows(mapply(function(gene_cluster) {
    gene_summary_stats <- extract_summary_statistics(gene_cluster, empirical_dataset = empirical_dataset, variant_reference = variant_reference)

    print(gene_summary_stats)
    message(sprintf("Input has %s variants", nrow(gene_summary_stats)))

    sample_size_threshold <- max(gene_summary_stats$sample_size) * min_sample_size_prop
    gene_summary_stats <- gene_summary_stats %>%
      filter(sample_size >= sample_size_threshold, i_squared <= max_i_squared)

    message(sprintf("After filtering %s variants remaining", nrow(gene_summary_stats)))

    if (normalize_sumstats) {
      message("Normalising effect sizes...") 

      min_sample_size <- min(gene_summary_stats %>% pull(sample_size))

      gene_summary_stats <- gene_summary_stats %>% mutate(
        Z = beta / standard_error,
        beta = Z / sqrt(sample_size + Z^2),
        standard_error = 1 / sqrt(min_sample_size + Z^2),
        Z = beta / standard_error)
    }

    variant_order_filtered <- variant_order[variant_order %in% gene_summary_stats$variant_index]
    
    gene_summary_stats <- gene_summary_stats[match(variant_order_filtered, (gene_summary_stats$variant_index)), ]
    gene_summary_stats <- gene_summary_stats[gene_summary_stats$variant_index %in% variant_order_filtered, ]
    print(gene_summary_stats)
    # Do fine-mapping
    if(nrow(gene_summary_stats) > 0 & all(gene_summary_stats$variant_index == variant_order_filtered) & !dry_run){

      estimated_res_var <- T
      fitted_rss2 <- tryCatch({
        fitted_rss2 <- susie_rss(bhat = gene_summary_stats$beta, shat = gene_summary_stats$standard_error, R = as.matrix(ld_matrix[variant_order_filtered, variant_order_filtered]), n = max(gene_summary_stats$sample_size), L = nCS, estimate_residual_variance = T, verbose=T)
      }, error = function(e) {
        fitted_rss2 <- susie_rss(bhat = gene_summary_stats$beta, shat = gene_summary_stats$standard_error, R = as.matrix(ld_matrix[variant_order_filtered, variant_order_filtered]), n = max(gene_summary_stats$sample_size), L = nCS, estimate_residual_variance = F, verbose=T)
        estimated_res_var <- F
        return(fitted_rss2)
      })

      if(!fitted_rss2$converged){
        fitted_rss2 <- susie_rss(bhat = gene_summary_stats$beta, shat = gene_summary_stats$standard_error, R = as.matrix(ld_matrix[variant_order_filtered, variant_order_filtered]), n = max(gene_summary_stats$sample_size), L = nCS, estimate_residual_variance = F, verbose=T)
        estimated_res_var = F
      }

      print("Finished!")
      print(fitted_rss2$converged)

      if(fitted_rss2$converged){
        print(summary(fitted_rss2))
        gene_summary_stats$SusieRss_lambda = estimate_s_rss(z=gene_summary_stats$beta / gene_summary_stats$standard_error, R = as.matrix(ld_matrix[variant_order_filtered, variant_order_filtered]),n=max(gene_summary_stats$sample_size))
        gene_summary_stats$SusieRss_pip = fitted_rss2$pip
        gene_summary_stats$SusieRss_CS = NA
        gene_summary_stats$SusieRss_ResVar = estimated_res_var
        if(length(fitted_rss2$sets[[1]])!=0){
          for(l in 1:length(fitted_rss2$sets[[1]])){
            gene_summary_stats$SusieRss_CS[fitted_rss2$sets[[1]][[l]]]=gsub("L","",names(fitted_rss2$sets[[1]])[l])
          }
        }
        lbfOut = t(fitted_rss2$lbf_variable)
        if(ncol(lbfOut)<nCS){
          lbfOut = as.data.frame(lbfOut)
          for(j in 1:(nCS - ncol(lbfOut))){
            lbfOut[paste("lbf_cs",j,sep="_")] = NA
          }
        }
        colnames(lbfOut) = paste("lbf_cs",1:nCS,sep="_")
        gene_summary_stats = cbind(gene_summary_stats, lbfOut)
      } else {
        print("Did not converge")
        gene_summary_stats$SusieRss_pip = NA
        gene_summary_stats$SusieRss_CS = NA
        gene_summary_stats$SusieRss_ResVar = NA
        for(j in 1:nCS){
          gene_summary_stats[[paste("lbf_cs",j,sep="_")]] = NA
        }
      }

    } else {
      if (!dry_run) {
        print("Was not able to start fine-mapping: variant_order does not align with gene_summary_stats")
      } else {
        print("Dry-run enabled. Skipped Fine-mapping.")
      }
      gene_summary_stats$SusieRss_pip = NA
      gene_summary_stats$SusieRss_CS = NA
      gene_summary_stats$SusieRss_ResVar = NA
      for(j in 1:nCS){
        gene_summary_stats[[paste("lbf_cs",j,sep="_")]] = NA
      }
    }

    if (debug != FALSE & debug == 'RSparsePro') {
      print("Writing RSparsePro tables")
      debug_table <- gene_summary_stats %>%
        mutate(Z = beta / standard_error, P=2*pnorm(q=abs(Z), lower.tail=FALSE)) %>%
        rename(RSID = variant_index) %>%
        inner_join(variant_reference, by = c("RSID"="variant_index"))
      fwrite(debug_table, sprintf("summary_stats_%s.txt", gene), sep="\t", col.names=T, row.names=F, quote=F)
      fwrite(debug_table %>% select(RSID, Z), sprintf("Z_%s_RSparsePro.txt", gene), sep="\t", col.names=T, row.names=F, quote=F)
      fwrite(as.matrix(ld_matrix[variant_order_filtered, variant_order_filtered]), sprintf("LD_%s_RSparsePro.txt", gene), sep="\t", col.names=F, row.names=F, quote=F)
    }
    print(gene_summary_stats)
    return(gene_summary_stats)
  }, locus_bed %>% group_by(gene, gene_cluster) %>% group_split(), SIMPLIFY =F))

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
  cat("LD data directory:", args$ld, "\n")
  cat("LD type:", args$ld_type, "\n")
  cat("Empirical data directory:", args$empirical, "\n")
  cat("Variant reference file:", args$variant_reference, "\n")
  cat("Uncorrelated genes file:", args$uncorrelated_genes, "\n")
  cat("BED files:", paste(args$bed_files, collapse = ", "), "\n")
  cat("Max I2:", paste(args$max_i2, collapse = ", "), "\n")
  cat("Min sample size prop:", paste(args$min_n_prop, collapse = ", "), "\n")
  cat("No adjust stats:", paste(args$no_adjust_stats, collapse = ", "), "\n")
  cat("Debugging:", args$debug, "\n")
  cat("Dry-run:", args$dry_run, "\n")

  # get path of parquet from arguments
  sumstats_path <- args$empirical
  permuted_dataset_path <- args$ld
  ld_type <- args$ld_type
  variant_reference_path <- args$variant_reference

  variant_reference <- arrow::read_parquet(variant_reference_path)
  uncorrelated_genes <- fread(args$uncorrelated_genes, header=F)$V1

  normalize_sumstats <- !args$no_adjust_stats
  min_sample_size_prop <- args$min_n_prop
  max_i_squared <- args$max_i2

  dry_run <- args$dry_run
  debug <- args$debug

  # load in parquet with Robert's strategy
  # Assume parquet dataset
  ds_format = ifelse(str_ends(sumstats_path, fixed(".feather")), "feather", "parquet")
  empirical_dataset <- open_dataset(sumstats_path, format=ds_format)
  print(empirical_dataset)

  # Switched to new ld reference dataset files with just phenotypes, variant_indices, and rho values
  # Added filtering on uncorrelated genes, since the dataset is no longer prefiltered on uncorrelated genes 
  # (the entire set of permuted genes is in this dataset)
  if (ld_type == 'gene-set' && !is.null(args$uncorrelated_genes)) {
    uncorrelated_genes <- fread(args$uncorrelated_genes, header=F)$V1
    permuted_dataset <- arrow::open_dataset(permuted_dataset_path) %>% filter(phenotype %in% uncorrelated_genes)
    ld_func <- function(variant_index_start, variant_index_end) {
      return(get_ld_matrix_alt(permuted_dataset, variant_index_start, variant_index_end))
    }
  } else if (ld_type == 'pcs') {
    permuted_dataset <- arrow::open_dataset(permuted_dataset_path)
    ld_func <- function(variant_index_start, variant_index_end) {
      return(get_ld_matrix_wide(permuted_dataset, variant_index_start, variant_index_end))
    } 
  } else if (ld_type == 'dosages') {
    dosages <- arrow::open_dataset(permuted_dataset_path)
    message("Opened dosage LD parquet")
    ld_func <- function(variant_index_start, variant_index_end) {
      message(sprintf("Calculating LD for %s-%s", variant_index_start, variant_index_end))
      dosages_mat <- dosages %>%
        filter(between(variant_index, variant_index_start, variant_index_end)) %>%
        collect() %>% as.data.table() %>% as.matrix(rownames='variant_index')
      message("Collecting matrix done!")
      message(dim(dosages_mat))
      str(dosages_mat)

      dosages_mat <- dosages_mat - rowMeans(dosages_mat)
      # Standardize each variable
      dosages_mat <- dosages_mat / sqrt(rowSums(dosages_mat^2))
      # Calculate correlations
      ld_matrix <- tcrossprod(dosages_mat)
      return(ld_matrix)
    }
  } else if (ld_type == 'feather') {
    ld_func <- function(variant_index_start, variant_index_end) {
      message(sprintf("Reading LD for %s-%s", variant_index_start, variant_index_end))
      ld_file <- file.path(
        permuted_dataset_path,
        sprintf("LD_%s-%s.feather", variant_index_start, variant_index_end))
      message(sprintf("Reading from %s", ld_file))
      ld_matrix <- arrow::read_feather(ld_file) %>%
        as.data.table() %>%
        as.matrix(rownames='__index_level_0__')
      return(ld_matrix)
    }
  }

  fine_mapping_results_per_locus <- mapply(function(bed_file) {
    locus_bed <- fread(bed_file, col.names = c("chromosome", "start", "end", "gene", "gene_cluster", "cluster"))
    message(sprintf("Starting finemapping in %s:%s-%s (%s genes)", locus_bed$chromosome[1], min(locus_bed$start), max(locus_bed$end), nrow(locus_bed)))
    # Assumes fine_mapping_output is some sort of data.frame/data.table or tibble
    fine_mapping_output <- finemap_locus(empirical_dataset=empirical_dataset,
                                         ld_func=ld_func,
                                         locus_bed=locus_bed,
                                         variant_reference=variant_reference,
                                         min_sample_size_prop=min_sample_size_prop,
                                         max_i_squared=max_i_squared,
                                         normalize_sumstats=normalize_sumstats,
                                         debug=debug,
                                         dry_run=dry_run)
    return(fine_mapping_output)
  }, args$bed_files, SIMPLIFY = F)

  combined_results <- bind_rows(fine_mapping_results_per_locus)

  if (!dry_run) {
    fwrite(combined_results, "finemapped.results.tsv", sep="\t", quote=F, row.names=F, col.names=T)
  }
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
