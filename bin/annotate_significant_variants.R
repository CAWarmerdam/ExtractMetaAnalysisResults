#!/usr/bin/env Rscript

library(argparse)
library(tidyverse)
library(arrow)
library(data.table)


# Main function to annotate significant variants
annotate_variants <- function(signSubset, variantReference, outputPrefix) {
  # Load significant variants using fread from data.table (faster for large CSVs)
  sign_variants <- fread(signSubset)
  
  # Load variant reference using arrow and read_parquet
  variant_ref <- read_parquet(variantReference)
  
  # Merge the significant variants with the variant reference using 'variant_index' as the key
  merged_data <- sign_variants %>%
    inner_join(variant_ref, by = "variant_index") %>%
    as.data.table()  # Convert back to a regular dt after lazy operations

  chromosomes <- unique(merged_data$chromosome)
  
  # Write the merged data to a CSV file
  for (cur_chromosome in chromosomes) {
      outputFile <- sprintf("%s.chr%d.csv", outputPrefix, cur_chromosome)
      fwrite(merged_data %>% filter(chromosome == cur_chromosome), outputFile, col.names=T, row.names=F, quote=F, sep="\t")
      message("Merged data has been successfully written to ", outputFile)
  }
}


# Main script execution
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = TRUE)
  }

  parser <- ArgumentParser(description = "Annotate significant variants by merging with a reference dataset.")
  
  parser$add_argument("--sign-subset", type = "character",
                      help = "Path to the gzipped CSV file with significant variants")
  parser$add_argument("--variant-reference", type = "character",
                      help = "Path to the parquet file containing variant reference data")
  parser$add_argument("--output-prefix", type = "character",
                      help = "Path to the output CSV file where merged data will be written")
  
  args <- parser$parse_args(argv)

  annotate_variants(args$sign_subset, args$variant_reference, args$output_prefix)
}

# Run the script
main()

