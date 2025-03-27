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
path_dirs <- strsplit(Sys.getenv("PATH"), .Platform$path.sep)[[1]]
script_name <- "run_susie_over_loci.R"

for (dir in path_dirs) {
  possible_path <- file.path(dir, script_name)
  if (file.exists(possible_path)) {
    script_path <- possible_path
    break
  }
}

source(script_path)

library(CARMA)

# Declare constants

# Declare functions
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

  ld_matrix <- ld_func(variant_index_start, variant_index_end)
  message("Starting to calculate LD...")
  start.time <- Sys.time()
  ld_matrix <- ld_func(variant_index_start, variant_index_end)
  print(as_tibble(ld_matrix))
  variant_order <- rownames(ld_matrix)
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  message("LD calculation done!")
  print(time.taken)

  if (debug == 'locus-wise') {
    saveRDS(as.matrix(ld_matrix), sprintf("LD_%s:%d-%d.rds", locus_chromosome, locus_start, locus_end))
  }

  # Do the remaining bit for finemapping the locus
  fine_mapping_results <- bind_rows(mapply(function(gene, start, end) {

    # Get the variants to load
    gene_locus_variant_reference <- variant_reference %>%
      filter(chromosome == locus_chromosome, between(bp, start, end)) %>% collect()

    # Get the indices of the variants to laod
    gene_variant_index_start <- min(gene_locus_variant_reference$variant_index)
    gene_variant_index_end <- max(gene_locus_variant_reference$variant_index)

    gene_summary_stats <- empirical_dataset %>%
      filter(phenotype == gene, between(variant_index, gene_variant_index_start, gene_variant_index_end)) %>%
      as.data.table()
    
    sample_size_threshold <- max(gene_summary_stats$sample_size) * min_sample_size_prop
    gene_summary_stats <- gene_summary_stats %>%
      filter(sample_size >= sample_size_threshold, i_squared <= max_i_squared)

    if (normalize_sumstats) {
      gene_summary_stats <- gene_summary_stats %>% mutate(
        Z = beta / standard_error,
        beta = Z / sqrt(sample_size + Z^2),
        standard_error = 1 / sqrt(sample_size + Z^2)
      )
    }

    variant_order_filtered <- variant_order[variant_order %in% gene_summary_stats$variant_index]
    
    gene_summary_stats <- gene_summary_stats[match(variant_order_filtered, (gene_summary_stats$variant_index)), ]
    gene_summary_stats <- gene_summary_stats[gene_summary_stats$variant_index %in% variant_order_filtered, ]
    print(gene_summary_stats)
    # Do fine-mapping
    if(all(gene_summary_stats$variant_index == variant_order_filtered) & !dry_run){

      z.list<-list()
      z.list[[1]]<-(gene_summary_stats$beta / gene_summary_stats$standard_error)
      ld.list<-list()
      ld.list[[1]]<-as.matrix(ld_matrix[variant_order_filtered, variant_order_filtered])
      lambda.list<-list()
      lambda.list[[1]]<-1
      CARMA.result<-CARMA(z.list,ld.list=ld.list,lambda.list = lambda.list, outlier.switch=T)

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
  }, locus_bed$gene, locus_bed$start, locus_bed$end, SIMPLIFY =F))

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
  empirical_dataset <- open_dataset(sumstats_path)

  # Switched to new ld reference dataset files with just phenotypes, variant_indices, and rho values
  # Added filtering on uncorrelated genes, since the dataset is no longer prefiltered on uncorrelated genes 
  # (the entire set of permuted genes is in this dataset)
  if (ld_type == 'gene-set') {
    permuted_dataset <- arrow::open_dataset(permuted_dataset_path) %>% filter(phenotype %in% uncorrelated_genes)
    ld_func <- function(variant_index_start, variant_index_end) {
      return(get_ld_matrix_alt(permuted_dataset, variant_index_start, variant_index_end))
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
  }

  fine_mapping_results_per_locus <- mapply(function(bed_file) {
    locus_bed <- fread(bed_file, col.names = c("chromosome", "start", "end", "gene", "cluster"))
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
