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
  default=0.5,
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

get_ld_matrix <- function(permuted_dataset, variant_indices) {
  rho_mat <- permuted_dataset %>%
    filter(variant_index %in% variant_indices) %>%
    collect() %>% as.data.table() %>%
    as.matrix(rownames=1)

  return(ld_calculator(rho_mat))
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

prune_susie_output <- function(chromosome_eqtl_table, ld_func, r2_threshold) {
  chromosome_eqtl_table <- chromosome_eqtl_table %>%
    group_by(variant_index) %>%
    slice_max(max_lbf)

  if (nrow(chromosome_eqtl_table) == 1) {
    return(chromosome_eqtl_table)
  }

  # Get the LD matrix
  ld_matrix <- ld_func(chromosome_eqtl_table$chromosome[1], chromosome_eqtl_table$variant_index)
  ld_matrix <- ld_matrix[
    as.character(chromosome_eqtl_table$variant_index),
    as.character(chromosome_eqtl_table$variant_index)]

  if (!all(rownames(ld_matrix) == as.character(chromosome_eqtl_table$variant_index))) {
  stop("LD matrix row names must match variant IDs in variants table")
  }

  # Create new columns for tagging information
  chromosome_eqtl_table$tagging_variants <- NA
  chromosome_eqtl_table$tagging_r <- NA

  # Sort variants by max_lbf (descending)
  chromosome_eqtl_table <- chromosome_eqtl_table %>% arrange(desc(max_lbf))

  # Keep track of removed variants
  to_remove <- rep(FALSE, nrow(chromosome_eqtl_table))
  names(to_remove) <- rownames(ld_matrix)

  # Iterate through variants
  for (index_variant in as.character(chromosome_eqtl_table$variant_index)) {
    if (to_remove[index_variant]) next  # Skip already tagged variants

    ld_values <- ld_matrix[index_variant, ]  # Get r^2 values

    # Identify variants in high LD
    which_variants_high_ld <- which(
      (ld_values^2) > r2_threshold &
        names(ld_values) != index_variant)
    high_ld_variants <- ld_values[which_variants_high_ld]

    if (length(high_ld_variants) > 0) {
      # Tag these variants
      chromosome_eqtl_table$tagging_variants[chromosome_eqtl_table$variant_index == as.integer(index_variant)] <- paste0(names(high_ld_variants), collapse=",")
      chromosome_eqtl_table$tagging_r[chromosome_eqtl_table$variant_index == as.integer(index_variant)] <- paste0(unname(high_ld_variants), collapse=",")

      # Mark them for removal
      to_remove[names(high_ld_variants)] <- TRUE
    }
  }

  return(chromosome_eqtl_table)
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

  ld_path <- "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/output/ld_panel_10497_pcs/"
  ld_type <- "pcs"
  susie_out_path <- "independent_variants_filtered_lbf2_mlog10p5_annotated_20250509.txt.gz"
  lead_variants <- "lead_variants.txt.gz"
  lower_r2_threshold <- 0.5
  variant_reference_path <- "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet"
  crossmapping_filter <- "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/2025-04-14-transCrossMap/output_20250509/10mbcis-35bpreads-10bpshift.txt.gz"
  cis_window <- 1e6
  twilight_window <- 5e6

  # Access the arguments
  ld_path <- args$ld
  cat("LD data directory:", ld_path, "\n")
  ld_type <- args$ld_type
  cat("LD type:", ld_type, "\n")
  susie_out_path <- args$independent_variants
  cat("Input table:", susie_out_path, "\n")
  lower_r2_threshold <- args$r2_threshold
  cat("R2 threshold", lower_r2_threshold, "\n")
  cat("Variant reference":variant_reference_path)

  variant_reference_path <- args$variant_reference

  # gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
  #   filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand)

  # hla_genes <- gene_reference %>%
  #   filter(
  #     (seqid == 6 & start <= 34000000 & 25000000 <= end)
  #   ) %>%
  #   pull(gene_id)
  #
  # fread(lead_variants) %>%
  #   filter(phenotype %in% hla_genes) %>%


  # Load crossmapping table
  crossmapping_output <- fread(crossmapping_filter) %>%
    mutate(crossmapping_flag = Proportion > 0.05,
           variant = str_replace_all(Variant, "_", ":"),
           variant = case_when(str_detect(Variant, "HGSV") ~ Variant, TRUE ~ variant)) %>%
    select(c("variant", "phenotype" = "Gene", "crossmapping_flag"))

  # get path of parquet from arguments
  eqtl_table <- fread(susie_out_path) %>% left_join(crossmapping_output) %>%
    filter(is.na(crossmapping_flag) | !crossmapping_flag)

  variant_reference <- arrow::read_parquet(variant_reference_path)

  permuted_dataset_path <- ld_path
  r_squared_threshold <- lower_r2_threshold

  # Switched to new ld reference dataset files with just phenotypes, variant_indices, and rho values
  # Added filtering on uncorrelated genes, since the dataset is no longer prefiltered on uncorrelated genes
  # (the entire set of permuted genes is in this dataset)
  permuted_dataset <- arrow::open_dataset(permuted_dataset_path)
  ld_func <- function(chromosome, variant_indices) {
    return(get_ld_matrix(permuted_dataset %>% filter(chr==chromosome), variant_indices))
  }

  ld_func_window <- function(chromosome, variant_index_start, variant_index_stop) {
    rho_mat <- permuted_dataset %>% filter(
        chr == chromosome,
        between(variant_index, variant_index_start, variant_index_stop)) %>%
      collect() %>% as.data.table() %>%
      as.matrix(rownames=1)

    return(ld_calculator(rho_mat))
  }

  precalculated_ld_path <- "../ld_matrices_20250509.rds"
  if (file.exists(precalculated_ld_path)) {
    ld_matrices <- readRDS(precalculated_ld_path)
  } else {
    variants_per_chromosome <- eqtl_table %>%
      group_by(chromosome) %>%
      distinct(variant_index) %>%
      group_split()

    names(variants_per_chromosome) <- sapply(variants_per_chromosome, function(x) x$chromosome[1])

    ld_matrices <- mapply(function(chrom) {
      ld_func(as.integer(chrom), variants_per_chromosome[[chrom]]$variant_index)
    }, names(variants_per_chromosome), SIMPLIFY=F, USE.NAMES = T)

    saveRDS(ld_matrices, precalculated_ld_path)
  }

  ld_func_mat <- function(chromosome, variant_indices) {
    chromosome_ld_matrix <- ld_matrices[[as.character(chromosome)]]
    return(chromosome_ld_matrix[as.character(variant_indices), as.character(variant_indices)])
  }

  pruned <- bind_rows(mapply(function(.x) {
    prune_susie_output(.x, ld_func_mat, 0.9)
  }, eqtl_table %>%
    group_by(chromosome, phenotype) %>%
    group_split(), SIMPLIFY=F)) %>%
    mutate(across(c(tagging_variants, tagging_r)
      , ~ lapply(str_split(.x, ","), as.numeric))) %>%
    group_by(chromosome, phenotype) %>%
    mutate(high_r2_flag = variant_index %in% c(unlist(tagging_variants)))

  fwrite(pruned, "independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-crossmapping-maxR2_0.9.txt.gz",
         sep="\t", row.names=F, col.names=T, quote=T)
  # fwrite(pruned, "independent_variants_filtered_lbf2_mlog10p5_annotated_20250415_filtered-crossmapping-maxR2_0.5.txt.gz",
  #        sep="\t", row.names=F, col.names=T, quote=T)
  #
  # independent_signals <- bind_rows(mapply(function(.x) {
  #   prune_susie_output(.x, ld_func_mat, 0.1)
  # }, eqtl_table %>% group_by(chromosome) %>% group_split(), SIMPLIFY=F)) %>%
  #   mutate(across(c(tagging_variants, tagging_r)
  #     , ~ lapply(str_split(.x, ","), as.numeric))) %>%
  #   group_by(chromosome) %>%
  #   filter(!variant_index %in% c(unlist(tagging_variants)))

  # Load processed set of independent variants
  # - Remove broad HLA region
  # - Annotate with variant class
  independent_variant_set_annot <- pruned %>%
    mutate(
      hla_flag = (chromosome == 6 & between(bp,25000000, 34000000)),
      variant_class = case_when(
        str_starts(variant, "HGSV") ~ "Structural variant",
        str_length(eff_allele) == 1 & str_length(non_eff_allele) == 1 ~ "SNP",
        TRUE ~ "InDel"),
      type = case_when(
        same_chromosome & between(bp, start - cis_window, end + cis_window) ~ "cis",
        same_chromosome & between(bp, start - twilight_window, end + twilight_window) ~ "intermediate",
        TRUE ~ "trans"
      ))

  print(independent_variant_set_annot %>%
          group_by(hla_flag, high_r2_flag, type) %>%
          tally())

  independent_variant_set_annot_filtered <- independent_variant_set_annot %>% filter(!hla_flag, !high_r2_flag)

  fwrite(independent_variant_set_annot_filtered,
         "independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping.txt.gz",
         col.names=T, row.names=F, sep="\t", quote=F)

  fwrite(pruned01,
         "independent_variants_filtered_lbf2_mlog10p5_annotated_20250827_filtered-maxR2_0.1-noHla-noCrossmapping.txt.gz",
         col.names=T, row.names=F, sep="\t", quote=F)

  independent_variant_set_annot_filtered <- fread("independent_variants_filtered_lbf2_mlog10p5_annotated_20250415_filtered-maxR2_0.9-noHla-noCrossmapping.txt.gz")
  # Load all variants that ended up in a credible set
  value <- independent_variant_set_annot_filtered %>%
    select(variant_index, phenotype, cluster, gene_cluster, CS_identifier)

  # For all credible set variants, keep those that are retained after filtering on HLA, crossmapping and pruned on R2
  finemapping_per_variant <- fread(variables["finemapping_results"])

  finemapping_per_variant_proc <- finemapping_per_variant %>%
    select(variant_index, phenotype, cluster, gene_cluster, SusieRss_CS, SusieRss_pip) %>%
    left_join(value) %>%
    inner_join(variant_reference) %>%
    group_by(chromosome, cluster, gene_cluster, SusieRss_CS, phenotype) %>%
    mutate(credible_set_retained = any(!is.na(CS_identifier)),
           lead_variant_in_cs = !is.na(CS_identifier))

  cs_variants <- finemapping_per_variant_proc %>% group_by(chromosome, cluster, gene_cluster, SusieRss_CS, phenotype) %>%
    arrange(desc(lead_variant_in_cs)) %>%
    summarise(CS_variant_index_list = paste(variant_index, collapse=","), CS_variant_list = paste(variant, collapse=","), variant_index_lead = variant_index[1])

  independent_variant_set_annot_filtered_with_cs_variants <- left_join(independent_variant_set_annot_filtered, cs_variants, by = c("variant_index" = "variant_index_lead", "chromosome", "cluster", "gene_cluster", "phenotype")) %>%
    mutate(CS_variant_list = ifelse(is.na(pip), as.character(variant), CS_variant_list),
           CS_variant_index_list = ifelse(is.na(pip), as.character(variant_index), CS_variant_index_list))

  fwrite(independent_variant_set_annot_filtered_with_cs_variants, "independent_variants_filtered_lbf2_mlog10p5_annotated_20250415_filtered-maxR2_0.9-noHla-noCrossmapping_with_CS_variants.txt.gz", sep="\t" ,row.names=F, col.names=T, quote=F)

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}
