#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(readxl)
library(data.table)
library(dtplyr)
library(arrow)
library(rtracklayer)
library(extrafont)
library(tidyr)
library(broom)
library(ggupset)

# Declare constants
loadfonts()
old <- theme_set(theme_classic(base_size = 10))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  text = element_text(family="Helvetica"),
  title = element_text(colour = "#595A5C", face = "bold"),
)

correction_vector <- c(
  "0 PCs per cohort (prime)", "0 PCs per cohort",
                       "20 PCs per cohort",
                       "25% of expl. var. per cohort",
                       "50% of expl. var. per cohort",
                       "60% of expl. var. per cohort",
                       "50 PCs per cohort")

# Assign names to the vector
names(correction_vector) <- c("empirical0", "empirical1", "empirical2", "empirical3", "empirical4", "empirical5", "empirical6")


zToBeta <- function(z, maf, n) {
  chi2 <- z^2
  a <- 2.0 * maf * (1.0 - maf) * (n + chi2)
  beta <- z / sqrt(a)
  se <- 1.0 / sqrt(a)
  return(list("beta" = beta, "se" = se))
}

# Declare function definitions
# gptRtoT <- function(r, n) (r*sqrt(n-2))/(sqrt(1-r^2))
# gptTtoR <- function(t, n) sqrt((t^2)/(n-2+t^2))
corToZ <- function(r, df){
  t <- sqrt(df) * abs(r) / sqrt(1 - (r*r))
  logp <- pmin(pt(t, df, log.p = T), pt(t, df, lower.tail=FALSE, log.p = T))#normaly need *2 since this is two-side. However only used to covert to z-score
  z <- qnorm(logp, log.p = T) #no /2 because because logp was not multipled by 2
  z <- ifelse(r > 0,z * -1,z)
  return(z)
}


zToR <- function(z,n) {
  r <-  z * 10.0^(0.00332727459306132 + -0.500803673176124 * log10(n)) - z^3.0 * 10.0^(-0.600314180613649 + -1.50039529106925 * log10(n)) + z^5.0 * 10.0^(-1.3415016045602 + -2.48579979106749 * log10(n)) - z^7.0 * 10.0^(-2.53051421275022 + -3.39996348904233 * log10(n))
  return(r)
}

zToCor <- function(z, df){
  t <- qt(pnorm(-abs(z),log.p = T), df, log.p=T)
  r <- t/sqrt(df+t^2)
  if(z > 0){
    r <- r * -1
  }
  return(r)
}

process_and_summarize <- function(data, filter_cell_type_snp = TRUE, subsample = FALSE) {
  max_finite_z_score <- max(data %>% ungroup() %>% filter(is.finite(z_score)) %>% pull(z_score))
  cat(max_finite_z_score)

  filtered_data <- data %>%
    group_by(set) %>% filter(!filter_cell_type_snp | !cell_type_snp) %>%
    mutate(z_score = if_else(is.infinite(z_score), max_finite_z_score, z_score))

  correlation_summary <- filtered_data %>% group_by(set, type) %>%
    summarise(tidy(cor.test(z_score, z_score_cell_type, use = "pairwise.complete.obs")))

  filtered_data <- filtered_data %>% filter(significant)

  n_eqtls_cis <- Inf
  n_eqtls_trans <- Inf
  n_subsample_table <- tibble(type = c("cis", "trans"), n_subsample = c(Inf, Inf))
  if (subsample) {
    n_eqtls_cis <- min(filtered_data %>% filter(type == "cis") %>% group_by(set) %>% summarise(n = n()) %>% pull(n))
    n_eqtls_trans <- min(filtered_data %>% filter(type == "trans") %>% group_by(set) %>% summarise(n = n()) %>% pull(n))
    n_subsample_table <- filtered_data %>% group_by(type, set) %>% summarise(n = n()) %>% group_by(type) %>% summarise(n_subsample = min(n))
  }
  print(n_eqtls_cis)
  print(n_eqtls_trans)
  print(n_subsample_table)

  processed_data <- filtered_data %>%
    inner_join(n_subsample_table, by ="type") %>%
    group_by(set, type) %>%
    filter(rank(-abs(z_score), ties.method = "random") <= n_subsample[1]) %>%
    mutate(
      n = n(),
      nom_sign = abs(z_score_cell_type) > -qnorm(0.05/2),
      concordant = sign(z_score) == sign(z_score_cell_type)
    )

  summarized_data <- processed_data %>% group_by(set, type) %>%
    summarise(
      tidy(cor.test(z_score, z_score_cell_type, use = "pairwise.complete.obs")),
              eqtls = sum(significant, na.rm=T),
              replicated = sum(nom_sign, na.rm=T),
              percent_replicated = replicated / eqtls * 100,
              concordant = sum(concordant & nom_sign, na.rm=T),
              percent_concordant = concordant / replicated * 100,
              cell_type_snps = n_distinct(variant[cell_type_snp]),
              percent_cell_type_snp_eqtls = sum(cell_type_snp, na.rm=T) / eqtls * 100,
              percent_cell_type_snps = cell_type_snps / n_distinct(variant) * 100) %>%
    inner_join(correlation_summary, by=c("set", "type"), suffix=c("_significant_discovery", "_all")) %>%
    mutate(message = sprintf(paste(
      "nr. of eQTLs: %s",
      "nominally replicated: %.0f (%.0f%%)",
      "concordant: %.0f (%.0f%%)",
      "cell type SNPs: %.0f (%.0f%%, %.0f%% of eQTLs)",
      "correlation of Z (confidence): %.2f (%.2f, %.2f)", sep="\n"),
      eqtls, replicated, percent_replicated, concordant, percent_concordant, cell_type_snps, percent_cell_type_snps, percent_cell_type_snp_eqtls, estimate_all, conf.low_all, conf.high_all))

  list(processed_data = processed_data, summarized_data = summarized_data)
}

create_ggplot <- function(data, cis_trans) {
  plot <- ggplot(data[[1]] %>% filter(nom_sign, type == cis_trans), aes(z_score, z_score_cell_type)) +
    geom_vline(aes(xintercept = 0), linewidth = (0.5 / (ggplot2::.pt * 72.27 / 96)), color = "black", alpha = 0.2) +
    geom_hline(aes(yintercept = 0), linewidth = (0.5 / (ggplot2::.pt * 72.27 / 96)), color = "black", alpha = 0.2) +
    geom_point(aes(color = concordant), shape = 16, size = 0.4, alpha = 0.4, show.legend = F) +
    scale_colour_manual(values = c("#D55E00", "black"), name = "Matching direction of effect") +
    geom_label(data = data[[2]] %>% filter(type == cis_trans),
               aes(label = message, x = -Inf, y = Inf),
               inherit.aes = FALSE,
               label.padding = unit(0.5, "lines"),
               label.r = unit(0, "lines"), label.size = 0, size = 6*0.35,
               vjust = "inward", hjust = "inward", colour = "black", fill = NA) +
    xlab("eQTLGen Z-score (meta-analyzed)") +
    ylab("sc-eQTLGen Z-score (meta-analyzed, oneK1K)") +
    facet_wrap(vars(set))

  return(plot)
}

get_replication_meta <- function(replication_path, replication_path_full) {
  if (file.exists(replication_path)) {
    replication_meta <- readRDS(replication_path)
  } else {
    replication <- fread(replication_path_full) %>%
      as_tibble() %>%
      pivot_longer(
        starts_with(c("p_value_", "beta_", "beta_Se_")),
        names_to = c(".value", "cell_type"),
        names_pattern = "(p_value|beta_Se|beta|n_samples|maf|hwe_p)_(.+)")

    replication_meta <- replication %>%
      group_by(EnsemblGene, SNP, assessed_allele) %>%
      mutate(weight = 1 / beta_Se^2,
             beta_weighted = beta * weight) %>%
      summarise(
        weight_summed = sum(weight, na.rm = T),
        se_meta_analysed = sqrt(1 / weight_summed),
        beta_meta_analysed = sum(beta_weighted, na.rm = T) / weight_summed) %>%
      mutate(
        z_score_cell_type = beta_meta_analysed / se_meta_analysed) %>%
      collect()

    rm(replication)

    saveRDS(replication_meta, replication_path)
  }
  return(replication_meta)
}

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  base_dir <- "/gpfs/space/GI/eQTLGen"
  base_dir <- "/Users/cawarmerdam/Documents/projects/eQTLGen"
  base_dir <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2"

  cell_type_composition_snps <- fread(file.path(base_dir, "public_data/gwas_catalog/processed/CellTypeCompSnpList.txt.gz"))

  pruned_snps <- fread(file.path(base_dir, "correction_optimization/pruned_immune_snps/immune_snps_chrAll.prune.in"),
                       header=F, col.names = c("variant")) %>%
    mutate(cell_type_snp = variant %in% cell_type_composition_snps$SNPName_eQTLGen) %>%
    collect()

  snp_ref <- read_parquet(file.path(base_dir, "processed_data/variants/1000G-30x_index.parquet")) %>%
    filter(variant %in% c(pruned_snps$variant))

  gene_ref <- file.path(base_dir, "public_data/Homo_sapiens.GRCh38.106.gtf.gz")
  gene_ref <- readGFF(gene_ref)
  gene_ref <- gene_ref[gene_ref$type == "gene", ]
  gene_ref <- unique(gene_ref[, c(1, 4, 5, 9, 11)])

  replication_path_full <- file.path(base_dir, "processed_data/sc-replication/20240830.combinedReplicationOneK1K_full.tsv.gz")
  replication_path <- file.path(base_dir, "processed_data/sc-replication/20240830.combinedReplicationOneK1K_full_meta.rds")

  replication_meta <- get_replication_meta(replication_path, replication_path_full)

  extraction_files <- Sys.glob(file.path(base_dir, "freeze3/Interpretation/extractions/correction_optimization/*/pruned/extracted_merged.txt.gz"))
  analysis_labels <- c("empirical0", "empirical1", "empirical2", "empirical3", "empirical4", "empirical5", "empirical6")
  named_strings <- analysis_labels %>%
    map_chr(~ extraction_files[grepl(.x, extraction_files)]) %>%
    set_names(analysis_labels)
  extractions <- bind_rows(mapply(fread, named_strings, SIMPLIFY=F, USE.NAMES = T), .id="set") %>%
    inner_join(snp_ref, by = c("variant_index"="variant_index")) %>%
    inner_join(gene_ref, by = c("phenotype" = "gene_id")) %>%
    mutate(chromosome = as.character(chromosome),
           type = case_when(chromosome == seqid & between(bp, start-5e6, end+5e6) ~ "cis",
                            TRUE ~ "trans"),
           significant = abs(z_score) > 5) %>%
    select(set, variant_index, variant, sample_size, c("gene_id" = "phenotype", "z_score", "beta", "significant", "non_eff_allele", "eff_allele", "type")) %>%
    mutate(set = correction_vector[set]) %>% collect()
  #
  # base_colnames <- c("variant_index", "gene_id", "R", "R_corrected", "sign", "sign_corrected")
  # colnames_vector <- c(
  #   paste0(base_colnames, ".100Comps"), "empty1",
  #   paste0(base_colnames, ".250Comps"), "empty2",
  #   paste0(base_colnames, ".500Comps")
  # )
  #
  # n_mat <- fread(file.path(base_dir, "freeze3/Interpretation/extractions/correction_optimization/empirical1_RNAseq_4GenPCTechCov_2024-06-26/exported_matrix/mat.n.txt.gz"), header=T) %>%
  #   as.matrix(rownames="variant")
  #
  # # Process input
  # global_correction <- read_excel(
  #   "/Users/cawarmerdam/Downloads/ImmuneSNPs-TransEQTLCorrectionResultsFor100-250-500PCs.xlsx",
  #   col_names=colnames_vector, skip=1)
  #
  # global_long <- global_correction %>%
  #   select(-empty1, -empty2) %>%
  #   pivot_longer(cols = everything(),
  #                names_to = c(".value", "set"),
  #                names_pattern = "(.+)\\.(.+)") %>%
  #   inner_join(snp_ref, by=c("variant_index" = "variant_index")) %>%
  #   inner_join(gtf2, by = c("gene_id" = "gene_id"))
  #
  # global_long_with_sign <- global_long %>%
  #   group_by(gene_id) %>%
  #   mutate(sample_size = n_mat[as.character(variant_index),gene_id[1]]) %>%
  #   ungroup() %>%
  #   mutate(z_score_uncorrected = corToZ(R, sample_size-1), z_score_corrected = corToZ(R_corrected, sample_size-1),
  #          significant_uncorrected = abs(z_score_uncorrected) > 5, significant_corrected = abs(z_score_corrected) > 5) %>%
  #   select(set, variant_index, variant, gene_id, z_score_uncorrected, z_score_corrected, significant_uncorrected, significant_corrected, non_eff_allele, eff_allele) %>%
  #   pivot_longer(cols=c(z_score_uncorrected, z_score_corrected, significant_uncorrected, significant_corrected), names_to = c(".value", "set2"),
  #                names_pattern = "(z_score|significant)_(.+)") %>%
  #   mutate(set = case_when(set2 == "uncorrected" ~ "0Comps",
  #                          TRUE ~ set)) %>%
  #   select(-set2) %>%
  #   distinct()

  results <- readRDS("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/eqtls_post_hoc_corrected_comparison_uncorrelated_genes_independent_snps_beta.rds")

  res_summarised_0 <- results %>%
    filter(abs(z_scores_corrected)>5) %>%
    mutate(variant_index = as.integer(as.character(variant_index))) %>%
    inner_join(snp_ref) %>%
    inner_join(pruned_snps) %>% collect()

  res_summarised <- res_summarised_0 %>%
    group_by(phenotype, variant_index, cis, cell_type_snp) %>%
    summarise(combination = list(sort(paste0(corrected_pcs, " Global PCs")), collapse = ",")) %>%
    mutate(variant_index = as.numeric(as.character(variant_index))) %>%
    group_by(combination, cis, cell_type_snp) %>% summarise(n = n()) %>%
    group_by(combination) %>%
    mutate(summed_n = sum(n)) %>%
    arrange(summed_n) %>%
    mutate(type = factor(ifelse(cis, "cis", "trans"), levels=c("cis","trans", ordered=T)))

  p <- ggplot(res_summarised, aes(x = combination, y=n, fill=cell_type_snp, alpha=type)) +
    geom_col() +
    scale_x_upset(scale_name="x", sets=paste0(c("20", "50", "100", "250", "500"), " Global PCs"), order_by="degree") +
    scale_fill_manual(values=c("#0072B2", "#D55E00")) +
    scale_alpha_manual(values=c(0.5, 1)) + facet_grid(rows = vars(type))
  ggsave("res_summarised.pdf", p, width=11, height=11)

  global_long_with_sign_uncorrected <- results %>%
    filter(corrected_pcs == 20) %>%
    dplyr::rename(gene_id = phenotype,
                  z_score = z_scores_uncorrected) %>%
    mutate(significant = abs(z_score) > 5,
           set = "0 global PCs",
           type = ifelse(cis, "cis", "trans"),
           variant_index = as.numeric(as.character(variant_index))) %>%
    inner_join(snp_ref) %>%
    select(set, variant_index, variant, gene_id, z_score, non_eff_allele, eff_allele, significant, type) %>% collect()

  global_long_with_sign_partial <- results %>%
    dplyr::rename(gene_id = phenotype,
                  z_score = z_scores_corrected_partial) %>%
    mutate(significant = abs(z_score) > 5,
           set = paste0(corrected_pcs, " global PCs"),
           type = ifelse(cis, "cis", "trans"),
           variant_index = as.numeric(as.character(variant_index))) %>%
    inner_join(snp_ref) %>%
    select(set, variant_index, variant, gene_id, z_score, non_eff_allele, eff_allele, significant, type) %>% collect()

  global_long_with_sign_corrected <- results %>%
    dplyr::rename(gene_id = phenotype,
           z_score = z_scores_corrected) %>%
    mutate(significant = abs(z_score) > 5,
           set = paste0(corrected_pcs, " global PCs (prime)"),
           type = ifelse(cis, "cis", "trans"),
           variant_index = as.numeric(as.character(variant_index))) %>%
    inner_join(snp_ref) %>%
    select(set, variant_index, variant, gene_id, z_score, non_eff_allele, eff_allele, significant, type) %>% collect()

  all_global <- bind_rows(global_long_with_sign_partial, global_long_with_sign_corrected, global_long_with_sign_uncorrected)

  summarised <- all_global %>% group_by(set) %>%
    summarise(n_pairs = n(), n_genes = n_distinct(gene_id), n_variants = n_distinct(variant), n_sign_cis = sum(significant[type == "cis"]), n_sign_trans = sum(significant[type == "trans"]))

  print(summarised, n =100)

  mj <- all %>% group_by(variant, gene_id) %>%
    filter(any(significant)) %>%
    distinct(variant, gene_id)

  write.table(mj, "correction_optimization_trans-eQTL_list.tsv.gz", sep="\t", row.names=F, col.names=T)

  gene_selection <- global_long_with_sign %>% nest_by(gene_id) %>% ungroup() %>% slice_sample(n=1000) %>% unnest(data) %>% distinct(gene_id)
  # Perform method
  # Process output
  fwrite(gene_selection, "gene_selection_immune_mediated_diseases.txt", col.names=T, row.names=F, sep="\t")

  # replication_monocytes <- replication %>%
  #   filter(cell_type == "CD4_T") %>%
  #   group_by(EnsemblGene, SNP, assessed_allele) %>%
  #   mutate(z_score_cell_type = beta / beta_se)

  bound_per_cohort <- as_tibble(extractions) %>% left_join(replication_meta, by = c("gene_id" = "EnsemblGene", "variant" = "SNP")) %>%
    mutate(z_score_cell_type = case_when(
      non_eff_allele == assessed_allele ~ -1*z_score_cell_type, TRUE ~ z_score_cell_type)) %>%
    inner_join(pruned_snps) %>% as_tibble()

  result_per_cohort <- process_and_summarize(bound_per_cohort, filter_cell_type_snp = F, subsample = F)
  saveRDS(result_per_cohort, "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/results_per_cohort.rds")

  bound_global <- as_tibble(all_global) %>% left_join(replication_meta, by = c("gene_id" = "EnsemblGene", "variant" = "SNP")) %>%
    mutate(z_score_cell_type = case_when(
      non_eff_allele == assessed_allele ~ -1*z_score_cell_type, TRUE ~ z_score_cell_type)) %>%
    inner_join(pruned_snps) %>% as_tibble()

  result_global <- process_and_summarize(bound_global, filter_cell_type_snp = F, subsample = F)
  saveRDS(result_global, "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/results_global_beta.rds")

  result_per_cohort <- readRDS("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/results_per_cohort.rds")

  # List of all lists
  all_lists <- list(result_global, result_per_cohort)

  # Combine all 'processed_data' tibbles
  combined_processed_data <- bind_rows(map(all_lists, "processed_data"))

  # Combine all 'summarized_data' tibbles
  combined_summarized_data <- bind_rows(map(all_lists, "summarized_data"))

  result2 <- process_and_summarize(combined_processed_data %>% filter(!str_ends(set, fixed("global PCs"))), filter_cell_type_snp = F, subsample = T)

  combined_wide <- combined_summarized_data %>%
    pivot_wider(set, names_from = "type", names_sep="_",
                values_from=c("eqtls", "replicated", "percent_replicated", "concordant", "percent_concordant",
                              "cell_type_snps", "percent_cell_type_snp_eqtls", "percent_cell_type_snps", "estimate", "estimate_significant_discovery"))

  fwrite(combined_wide, "comparison_summary_beta.csv", sep=",", row.names=F)

  result_all <- list(processed_data = combined_processed_data, summarized_data = combined_summarized_data)

  p1 <- create_ggplot(result_all, "trans")
  ggsave("correction_comparison_2024-09-02.pdf", p1, height=9, width=9)




  bound <- as_tibble(extractions) %>% left_join(replication_meta, by = c("gene_id" = "EnsemblGene", "variant" = "SNP")) %>%
    mutate(z_score_cell_type = case_when(
      non_eff_allele == assessed_allele ~ -1*z_score_cell_type, TRUE ~ z_score_cell_type)) %>%
    inner_join(pruned_snps) %>% as_tibble()

  result1 <- process_and_summarize(bound, filter_cell_type_snp = F, subsample = F)
  result2 <- process_and_summarize(bound, filter_cell_type_snp = F, subsample = T)
  result3 <- process_and_summarize(bound, filter_cell_type_snp = TRUE, subsample = T)

  p1 <- create_ggplot(result1)
  ggsave("correction_comparison_2024-08-20.pdf", p1, height=9, width=9)
  p2 <- create_ggplot(result2)
  ggsave( "correction_comparison_sampled_2024-08-06.pdf", p2, height=9, width=9)
  p3 <- create_ggplot(result3)
  ggsave( "correction_comparison_sampled_non_cell_type_2024-08-06.pdf", p3, height=9, width=9)

  bound_mono <- all %>% inner_join(replication_monocytes, by = c("gene_id" = "EnsemblGene", "variant" = "SNP")) %>%
    mutate(z_score_cell_type = case_when(
      non_eff_allele == assessed_allele ~ -1*z_score_cell_type, TRUE ~ z_score_cell_type)) %>%
    filter(significant) %>%
    mutate(set = factor(correction_vector[set], levels=unname(correction_vector)), ordered=T) %>%
    inner_join(pruned_snps)

  result1_mono <- process_and_summarize(bound_mono, filter_cell_type_snp = F, subsample = F)
  result2_mono <- process_and_summarize(bound_mono, filter_cell_type_snp = F, subsample = T)
  result3_mono <- process_and_summarize(bound_mono, filter_cell_type_snp = TRUE, subsample = T)

  p1_mono <- create_ggplot(result1_mono)
  ggsave("correction_comparison_mono_2024-08-06.pdf", p1_mono, height=9, width=9)
  p2_mono <- create_ggplot(result2_mono)
  ggsave( "correction_comparison_mono_sampled_2024-08-06.pdf", p2_mono, height=9, width=9)
  p3_mono <- create_ggplot(result3_mono)
  ggsave( "correction_comparison_mono_sampled_non_cell_type_2024-08-06.pdf", p3_mono, height=9, width=9)

  cell_type_comparison <- bind_rows(list(
    "with_cell_type_snps" = result2$summarized_data,
    "without_cell_type_snps"= result3$summarized_data), .id="analysis") %>%
    select(analysis, set, estimate, percent_replicated, percent_concordant) %>%
    pivot_wider(names_from="analysis", values_from=c("estimate", "percent_replicated", "percent_concordant"))

  fwrite(cell_type_comparison, "comparison_confinement_to_cell_type_snps.tsv", sep="\t")

  bound_wide <- bound %>%
    select(set, variant, gene_id, z_score, z_score_cell_type) %>%
    pivot_wider(id_cols = c(z_score_cell_type, variant, gene_id), names_from = "set", values_from = "z_score", names_prefix="z_score_")

  # Function to compute correlations
  compute_cor <- function(input_df, col1, col2, coly) {
    input_df <- na.omit(input_df %>% select(all_of(c(col1, col2, coly))))
    res <- cocor.dep.groups.overlap(
      cor.test(input_df[[col1]], input_df[[coly]])$estimate,
      cor.test(input_df[[col2]], input_df[[coly]])$estimate,
      cor.test(input_df[[col1]], input_df[[col2]])$estimate,
      n=nrow(input_df),
    test="hittner2003", return.htest = T)
    return(res$hittner2003$p.value)
  }

  col_pairs <- combn(colnames(bound_wide)[str_starts(colnames(bound_wide), "^z_score")], 2, simplify = FALSE)

  cor_results <- map_dfr(col_pairs, ~ {
    tibble(
      col1 = .x[1],
      col2 = .x[2],
      correlation = compute_cor(bound_wide, .x[1], .x[2], "z_score_cell_type")
    )
  })


  covariates <- fread("~/Documents/projects/eQTLGen/freeze3/eqtl_mapping/input/covariates_long_extended_2024-08-08.txt") %>%
    mutate(type = case_when(
      !is.na(PC) ~ "Expression PC",
      !is.na(str_match(covariate, "GenPC")) ~ "Genetic PC",
      !is.na(str_match(covariate, "intercept")) ~ "Intercept",
      TRUE ~ "Additional")) %>%
    group_by(cohort, type) %>%
    summarise(across(c(starts_with("analysis"), explained_variance), sum)) %>%
    ungroup() %>%
    arrange(analysis4_4GenPC50pExpPC) %>%
    mutate(cohort = ordered(cohort, levels=unique(cohort))) %>%
    pivot_longer(-c(cohort, type, explained_variance), names_to = "analysis", values_to="count") %>%
    group_by(analysis, type) %>%
    arrange(count) %>%
    mutate(index = row_number())


  ggplot(covariates, aes(count, index)) +
    geom_col(width=0.8, orientation='y') +
    facet_grid(rows = vars(analysis), cols=vars(type), scales="free_x", space="free_x") +
    theme(strip.text.x.top = element_text(angle = 90))
  ggsave("covariate_count_per_analysis.pdf", height=11, width=9)

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}