#!/usr/bin/env Rscript


# Load libraries
library(argparse)
library(data.table)
library(tidyverse)
library(arrow)

# Declare constants
argv <- c(
  "--input-parquet", "/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/output/all_empirical_per_cohort_4GenPC20ExpPC_2025-06-12/eqtls/cohort",
  "--gene-list", "ENSG00000272368", "ENSG00000272374", "ENSG00000272382",
  "--variant-reference", "/gpfs/space/GI/eQTLGen/processed_data/variants/1000G-30x_index.parquet",
  "--finemapped-eqtls", "/gpfs/space/GI/eQTLGen/fine_mapping/output/genome_wide_finemapping_susie_20250314/finemapped/fix_20250509/independent_variants_filtered_lbf2_mlog10p5_annotated_20250509_filtered-maxR2_0.9-noHla-noCrossmapping.txt.gz",
  "--maf-table", "/gpfs/space/GI/eQTLGen/freeze3/InclusionLists/inclusion_as_parquet/maf_table.parquet"
)


# Declare function definitions

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input

  # Create parser object
  parser <- ArgumentParser(description = 'Extract per-cohort results')

  # Define arguments
  parser$add_argument('--input-parquet', required = TRUE, help = 'Path to input Parquet file')
  parser$add_argument('--gene-list', nargs = '+', required = TRUE, help = 'List of genes')
  parser$add_argument('--maf-table', required = TRUE, help = 'Path to MAF table file')
  parser$add_argument('--variant-reference', required = TRUE, help = 'Path to variant reference file')
  parser$add_argument('--finemapped-eqtls', required = TRUE, help = 'Path to fine-mapped variants file')

  # Parse the arguments
  args <- parser$parse_args(argv)

  # Get gene list
  gene_list <- args$gene_list

  # Perform method
  # Load finemapped eQTLs
  finemapped_eqtls <- fread(args$finemapped_eqtls) %>%
    filter(phenotype %in% args$gene_list) %>%
    distinct(phenotype, variant, variant_index)

  # Get MAF
  maf_table <- open_dataset(args$maf_table) %>%
    filter(variant_index %in% finemapped_eqtls$variant_index) %>%
    collect() %>%
    pivot_longer(cols=-variant_index, names_to="cohort", values_to="cohort_eaf")

  # Load summary statistics and filter on finemapped eqtls, annotate with MAF, and annotate with sample size
  per_cohort_summary_statistics <- bind_rows(mapply(function(gene) {
    open_dataset(file.path(args$input_parquet, paste0("phenotype=", gene))) %>%
      collect() %>%
      filter(!is.na(beta)) %>%
      mutate(phenotype = gene) %>%
      inner_join(finemapped_eqtls, by = c("variant_index", "phenotype"))
  }, gene_list, SIMPLIFY=F)) %>%
    inner_join(maf_table, by =c("variant_index", "cohort"))

  # Process output
  filename_output <- sprintf(
    "extracted_finemapped_variants_per_cohort_with_gene_%s-and-%s-more.txt.gz",
    gene_list[1], length(gene_list) - 1)

  fwrite(per_cohort_summary_statistics, filename_output,
         sep="\t", row.names=F, col.names=T, quote=F)
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}