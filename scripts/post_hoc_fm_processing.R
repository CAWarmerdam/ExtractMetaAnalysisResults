#!/usr/bin/env Rscript


# Load libraries
library(data.table)
library(arrow)
library(tidyverse)
library(dtplyr)
library(dbplyr)
library(rtracklayer)
library(argparse)

parser <- ArgumentParser(
  description = "Process fine-mapping and reference files for eQTLGen phase 2 pipeline."
)

parser$add_argument(
  "--variant_reference", required=TRUE,
  help="Path to the variant reference parquet file."
)
parser$add_argument(
  "--gene_reference", required=TRUE,
  help="Path to the gene reference GTF (gzipped)."
)
parser$add_argument(
  "--finemapping_summary", required=TRUE,
  help="Path to the finemapping summary TSV file."
)
parser$add_argument(
  "--finemapping_results", required=TRUE,
  help="Path to the finemapping results TSV file."
)
parser$add_argument(
  "--gene_inclusion_file", required=TRUE,
  help="Path to the gene inclusion table (parquet)."
)
parser$add_argument(
  "--variant_inclusion_file", required=TRUE,
  help="Path to the variant inclusion table (parquet)."
)
parser$add_argument(
  "--maf_file", required=TRUE,
  help="Path to the MAF table (parquet)."
)
parser$add_argument(
  "--master_table", required=TRUE,
  help="Path to the master table (txt)."
)

# Declare constants
old <- theme_set(theme_classic(base_size = 10))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  text = element_text(family="Helvetica"),
  title = element_text(colour = "#595A5C", face = "plain")
)

# Declare function definitions

# Declare function definitions
get_maf_func <- function(gene_inclusion_file, variant_inclusion_file, maf_file, sample_sizes) {
  gene_inclusion_ds <- open_dataset(gene_inclusion_file)
  variant_inclusion_ds <- open_dataset(variant_inclusion_file)
  allele_frequency_ds <- open_dataset(maf_file)
  return(function(variant_indices, gene) {
    calculate_af(as.integer(variant_indices), gene, gene_inclusion_ds, variant_inclusion_ds, allele_frequency_ds, sample_sizes)
  })
}

calculate_af_per_gene <- function(variant_indices, gene, gene_inclusion_ds, variant_inclusion_ds, allele_frequency_ds, sample_sizes) {
  cohort_names <- names(sample_sizes)

  af <- setNames(rep(NA, length(variant_indices)), as.character(variant_indices))

  gene_inclusion_matrix <- gene_inclusion_ds %>% filter(phenotype == gene) %>%
    collect() %>% as.data.table() %>% as.matrix(rownames="phenotype")
  variant_inclusion_matrix <- variant_inclusion_ds %>% filter(variant_index %in% variant_indices) %>%
    collect() %>% as.data.table() %>% as.matrix(rownames="variant_index")
  af_matrix <- allele_frequency_ds %>% filter(variant_index %in% variant_indices) %>%
    collect() %>% as.data.table() %>% as.matrix(rownames="variant_index")

  variant_indices_as_character <- as.character(intersect(
    variant_indices, intersect(rownames(af_matrix), rownames(variant_inclusion_matrix))))

  gene_inclusion_vector <- gene_inclusion_matrix[gene, cohort_names, drop=F]
  af_matrix <- af_matrix[variant_indices_as_character, cohort_names, drop=F]
  variant_inclusion_matrix <- variant_inclusion_matrix[variant_indices_as_character, cohort_names, drop=F]

  af_matrix[!variant_inclusion_matrix] <- NA
  af_matrix[,!gene_inclusion_vector] <- NA

  presence_matrix <- variant_inclusion_matrix
  class(presence_matrix) <- "integer"

  af[variant_indices_as_character] <- colSums(t(af_matrix) * sample_sizes, na.rm=T) / colSums(t(presence_matrix) * sample_sizes)
  return(af)
}


calculate_af_high_mem <- function(variant_indices, genes, gene_inclusion_file, variant_inclusion_file, maf_file, sample_sizes) {
  cohort_names <- names(sample_sizes)

  gene_inclusion_ds <- open_dataset(gene_inclusion_file)
  variant_inclusion_ds <- open_dataset(variant_inclusion_file)
  allele_frequency_ds <- open_dataset(maf_file)

  gene_inclusion_matrix <- gene_inclusion_ds %>% filter(phenotype %in% genes) %>%
    collect() %>% as.data.table() %>% as.matrix(rownames="phenotype")
  variant_inclusion_matrix <- variant_inclusion_ds %>% filter(variant_index %in% variant_indices) %>%
    collect() %>% as.data.table() %>% as.matrix(rownames="variant_index")
  af_matrix <- allele_frequency_ds %>% filter(variant_index %in% variant_indices) %>%
    collect() %>% as.data.table() %>% as.matrix(rownames="variant_index")

  af_results <- vector("numeric", length(variant_indices))

  variant_indices_as_character <- as.character(variant_indices)

  # Loop over each predefined combination of variant and gene (both vectors are of equal length)
  for (i in seq_along(variant_indices)) {
    message(sprintf("%f/%f", i, length(variant_indices)))
    # Get the current variant and gene combination
    variant <- variant_indices_as_character[i]
    gene <- genes[i]

    # Calculate presence matrix for the current variant and gene
    presence <- gene_inclusion_matrix[gene, cohort_names, drop = TRUE] *
      variant_inclusion_matrix[variant, cohort_names, drop = TRUE]

    # Calculate the weighted allele frequency
    numerator <- presence * af_matrix[variant, cohort_names, drop = TRUE] * sample_sizes
    denominator <- presence * sample_sizes

    # Calculate final allele frequency, avoiding division by zero
    af_results[i] <- sum(numerator, na.rm = TRUE) / sum(denominator, na.rm = TRUE)
  }

  # Return results as a data.frame for easy viewing
  result_df <- data.frame(
    variant_index = variant_indices,
    gene = genes,
    weighted_AF = af_results
  )

  return(result_df)
}

variables <- function() {

  variables <- c(
    "variant_reference" = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet",
    "gene_reference" = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/public_data/Homo_sapiens.GRCh38.106.gtf.gz",
    "finemapping_summary" = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/output/genome_wide_finemapping_susie_20250314/finemapped/finemapping_summary.tsv",
    "finemapping_results" = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/fine_mapping/output/genome_wide_finemapping_susie_20250314/finemapped/finemapped.results.tsv",
    "gene_inclusion_file" = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/InclusionLists/inclusion_as_parquet/gene_inclusion_table.parquet",
    "variant_inclusion_file" = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/InclusionLists/inclusion_as_parquet/variant_inclusion_table.parquet",
    "maf_file" = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/InclusionLists/inclusion_as_parquet/maf_table.parquet",
    "master_table" = "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/eqtl_mapping/input/mastertable_empirical_2024-08-29_extended.txt"
  )

  variables <- c(
    "variant_reference" = "~/Documents/projects/eQTLGen/processed_data/variants/1000G-30x_index.parquet",
    "gene_reference" = "~/Documents/projects/eQTLGen/public_data/Homo_sapiens.GRCh38.106.gtf.gz",
    "finemapping_summary" = "~/Documents/projects/eQTLGen/fine_mapping/output/genome_wide_finemapping_susie_20250314/finemapped/finemapping_summary.tsv",
    "finemapping_results" = "~/Documents/projects/eQTLGen/fine_mapping/output/genome_wide_finemapping_susie_20250314/finemapped/finemapped.results.tsv",
    "gene_inclusion_file" = "~/Documents/projects/eQTLGen/freeze3/InclusionLists/inclusion_as_parquet/gene_inclusion_table.parquet",
    "variant_inclusion_file" = "~/Documents/projects/eQTLGen/freeze3/InclusionLists/inclusion_as_parquet/variant_inclusion_table.parquet",
    "maf_file" = "~/Documents/projects/eQTLGen/freeze3/InclusionLists/inclusion_as_parquet/maf_table.parquet",
    "master_table" = "~/Documents/projects/eQTLGen/freeze3/eqtl_mapping/input/mastertable_empirical_2024-08-29_extended.txt"
  )

}


# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  args <- parser$parse_args()

  # turn into a named character vector, just like your original variables <- c(...)
  variables <- unlist(args, use.names = TRUE)

  # Load sample sizes
  sample_sizes <- deframe(fread(variables["master_table"]) %>% select(cohort_new_name, dataset_n))
  # Get function to load Minor allele  frequencies
  maf_func <- get_maf_func(variables["gene_inclusion_file"], variables["variant_inclusion_file"], variables["maf_file"], sample_sizes)

  # Get variant reference
  variant_reference <- read_parquet(variables["variant_reference"])
  # Load gene reference
  gene_reference <- rtracklayer::readGFF(variables["gene_reference"]) %>%
    filter(type == "gene") %>% select(gene_id, gene_name, seqid, start, end, strand)

  # Read data using data.table for speed
  finemapping_summary <- fread(variables["finemapping_summary"]) %>%
    mutate(finemapping_index = row_number())

  # Load finemapping per variant (all finemapped variants that got assigned to a CS)
  finemapping_per_variant <- fread(variables["finemapping_results"]) %>%
    select(variant_index, beta, standard_error, phenotype, i_squared, sample_size, SusieRss_CS)

  # Windows for cis/trans/twilight eQTLs
  # Process input
  cis_window <- 1e6
  twilight_window <- 5e6

  setDT(variant_reference)

  variant_mapping_start <- variant_reference[
    , .(window_variant_index_start = min(variant_index)), by = .(chromosome, bp)
  ]

  variant_mapping_end <- variant_reference[
    , .(window_variant_index_end = min(variant_index)), by = .(chromosome, bp)
  ]

  collected <- fread("finemapping.log.collected.bed")

  all_windows <- collected %>% select(chromosome, gene, gene_cluster, cluster, start, end)
  all_windows <- collected %>% group_by(chromosome, gene, gene_cluster, cluster) %>% summarise(start = min(start), end = max(end))

  # Finemapping windows that failed
  finemapping_windows_failed <- finemapping_summary %>%
    filter(!converged | n_unique_CSs_pass == 0) %>%
    inner_join(variant_reference %>% select(chromosome, variant_index), by = c("min_variant_index" = "variant_index")) %>%
    inner_join(all_windows, by = c("chromosome", "phenotype" = "gene", "cluster", "gene_cluster")) %>%
    inner_join(variant_mapping_start, by = join_by(chromosome, closest(start <= bp))) %>%
    inner_join(variant_mapping_end, by = join_by(chromosome, closest(end >= bp)))

  meta_analysis <- open_dataset("/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/eqtl_mapping/output/2024-09-30_meta_analysis")

  # Per failed gene-window combination, get the lead variant.
  finemapping_windows_failed_queried <- pmap(
    finemapping_windows_failed,
    function(phenotype, window_variant_index_start, window_variant_index_end, finemapping_index, converged, ...) {
      print(sprintf("%s, %s: %s-%s", finemapping_index, phenotype, window_variant_index_start, window_variant_index_end))
      meta_analysis %>%
        filter(
          phenotype == !!phenotype,
          between(variant_index, !!window_variant_index_start, !!window_variant_index_end)
        ) %>%
        slice_max(abs(beta/standard_error), with_ties = F, n=1) %>% collect() %>%
        mutate(finemapping_index = finemapping_index, converged=converged)
    }
  )

  annotated_finemapping_failed <- bind_rows(finemapping_windows_failed_queried) %>%
    inner_join(finemapping_windows_failed, by = c(finemapping_index, phenotype, converge)) %>%
    select(-trace, -CS_size_pass, -lead_variant_standard_error, -lead_variant_index, -lead_variant_beta, -lead_variant_z, -window_variant_index_start, -window_variant_index_end) %>%
    inner_join(variant_reference, by = c("variant_index" = "variant_index")) %>%
    inner_join(gene_reference, by = c("phenotype" = "gene_id")) %>%
    mutate(same_chromosome = seqid == chromosome,
           type = case_when(
             same_chromosome & between(bp, start - cis_window, end + cis_window) ~ "cis",
             TRUE ~ "trans"
           )) %>%
    mutate(CS_identifier = sprintf("%s_NA", finemapping_index))

  renaming_vector <- c("variant_index" = "variant_indices", "p_value" = "p_values", "max_lbf" = "max_lbfs", "pip" = "pips")

  annotated_finemapping_succeeded <- finemapping_summary %>%
    filter(converged) %>%
    mutate(finemapping_index = row_number()) %>%
    mutate(
      across(c(variant_indices, pips, p_values, max_lbfs, CS_size), ~ lapply(str_split(.x, ","), as.numeric))) %>%
    unnest_longer(col = c(variant_indices, pips, p_values, max_lbfs, CS_size)) %>%
    rename(renaming_vector) %>%
    inner_join(variant_reference, by = c("variant_index" = "variant_index")) %>%
    inner_join(gene_reference, by = c("phenotype" = "gene_id")) %>%
    mutate(same_chromosome = seqid == chromosome,
           type = case_when(
             same_chromosome & between(bp, start - cis_window, end + cis_window) ~ "cis",
             TRUE ~ "trans"
           )) %>%
    inner_join(finemapping_per_variant,
               by = c("variant_index"="variant_index", "phenotype"="phenotype")) %>%
    mutate(CS_identifier = sprintf("%s_%s", finemapping_index, SusieRss_CS)) %>%
    select(-trace, -CS_size_pass, -lead_variant_standard_error, -lead_variant_index, -lead_variant_beta, -lead_variant_z)

  fwrite(annotated_finemapping_succeeded, "succeeded_finemapping_output_unfiltered_20250509.txt.gz", col.names=T, row.names=F, quote=F, sep="\t")

  borzoi_set <- annotated_finemapping_succeeded %>%
    filter(pip > 0.9, type == "cis", str_length(non_eff_allele) == 1, str_length(eff_allele) == 1)

  fwrite(borzoi_set, "borzoi_set_typeCis_pip0.9_snpsOnly.txt.gz", col.names=T, row.names=F, quote=F, sep="\t")

  independent_variant_set_unfiltered <- bind_rows(annotated_finemapping_failed, annotated_finemapping_succeeded)
  fwrite(annotated_finemapping_succeeded, "independent_variants_unfiltered_20250509.txt.gz", col.names=T, row.names=F, quote=F, sep="\t")
  independent_variant_set <- independent_variant_set_unfiltered %>%
    filter(max_lbf > 2 & p_value < 1e-5)

  af_tib <- calculate_af_high_mem(independent_variant_set$variant_index,
                                  independent_variant_set$phenotype,
                                  variables["gene_inclusion_file"],
                                  variables["variant_inclusion_file"],
                                  variables["maf_file"], sample_sizes)

  fwrite(af_tib, "allele_frequencies_independent_variants_all_20250509.txt.gz", col.names=T, row.names=F, quote=F, sep="\t")
  af_tib <- fread("allele_frequencies_independent_variants_all_20250509.txt.gz")

  independent_variant_set_annot <- independent_variant_set %>% inner_join(af_tib, by = c("variant_index", "phenotype"="gene")) %>%
    rename(all_of(c("eaf"="weighted_AF")))

  #fwrite(independent_variant_set_annot, "independent_variants_all_annotated_20250320.txt.gz", col.names=T, row.names=F, quote=F, sep="\t")

  fwrite(independent_variant_set_annot, "independent_variants_filtered_lbf2_mlog10p5_annotated_20250509.txt.gz", col.names=T, row.names=F, quote=F, sep="\t")

  unique_independent_variants <- independent_variant_set_annot %>%
    distinct(variant)
  fwrite(unique_independent_variants, "unique_independent_variants_filtered_lbf2_mlog10p5_annotated_20250509.txt.gz", col.names=F, row.names=F, quote=F, sep="\t")

  # Write table that is going to be the input to the crossmapping filter
  all_finemapped_variants <- independent_variant_set_annot %>%
    mutate(type2 = case_when(
             same_chromosome & between(bp, start - cis_window, end + cis_window) ~ "cis",
             same_chromosome & between(bp, start - twilight_window, end + twilight_window) ~ "twilight",
             TRUE ~ "trans"
           ),
           z_score = beta/standard_error
    ) %>% filter(type2 == "trans") %>%
    select(c("gene_id" = "phenotype"), variant, chromosome, bp, eff_allele, non_eff_allele, beta, standard_error, z_score)

  fwrite(all_finemapped_variants,
         "independent_variants_filtered_lbf2_mlog10p5_gene_variant_pairs_only_trans_5e6fromTssTes_20250509.txt.gz",
  sep="\t", col.names=T, row.names=F, quote=F)

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}