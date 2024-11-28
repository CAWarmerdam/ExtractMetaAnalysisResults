#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)

# Declare constants
loadfonts()
old <- theme_set(theme_classic(base_size = 10))
theme_update(
  line = element_line(
    colour = "black", linewidth = (0.5 / (ggplot2::.pt * 72.27/96)),
    linetype = 1, lineend = "butt", arrow = F, inherit.blank = T),
  strip.background = element_rect(colour = NA, fill = NA),
  panel.grid.minor = element_blank(),
  text = element_text(family="Helvetica"),
  title = element_text(colour = "#595A5C", face = "bold")
)


# Create a parser object
parser <- ArgumentParser(description = 'Correct trans-eQTLs for cis-eQTLs.')

# Add command-line arguments
parser$add_argument(
  "--input-file",
  required = TRUE,
  help = "Path to the directory containing empirical data."
)

parser$add_argument(
  "--genes",
  required = TRUE,
  help = "Genes to use"
)

parser$add_argument(
  "--cis-explained-variance",
  required = TRUE,
  help = "Genes to use"
)

parser$add_argument(
  "--variant-reference",
  required = TRUE,
  help = "Path to the variant reference file."
)

parser$add_argument(
  "--gene-reference",
  required = TRUE,
  help = "Path to the gene reference file."
)

parser$add_argument(
  "--output-prefix",
  required = TRUE,
  help = "Output prefix"
)



IdentifyLeadSNPs <- function(data, window = 1000000, Pthresh = 5e-8, snp_id_col = "SNP", snp_chr_col = "Chr", snp_pos_col = "Pos", beta_col = "beta", se_col = "se", p_col = NULL) {

  if (is.null(p_col)){
    data <- data.table(SNP = data[[snp_id_col]],
                       snp_chr = data[[snp_chr_col]],
                       snp_pos = data[[snp_pos_col]],
                       beta = data[[beta_col]],
                       se = data[[se_col]])

    data$P <- ZtoP(data$beta/data$se)
  } else {

    data <- data.table(SNP = data[[snp_id_col]],
                       snp_chr = data[[snp_chr_col]],
                       snp_pos = data[[snp_pos_col]],
                       beta = data[[beta_col]],
                       se = data[[se_col]],
                       P = data[[p_col]])

  }
  data$Z <- data$beta/data$se

  data_f <- data[data$P < as.numeric(Pthresh), ]
  # Iteratively identify most significant SNP, and remove all other SNPs in the window
  res <- data_f[-c(1:nrow(data_f)), ]

  while (min(data_f$P) <= Pthresh) {
    lead_snp <- data_f[abs(data_f$Z) == max(abs(data_f$Z)), ]
    if (nrow(lead_snp) > 1) {
      lead_snp <- lead_snp[1, ]
    }
    res <- rbind(res, lead_snp)
    data_f <- data_f[!(data_f$snp_chr == lead_snp$snp_chr & data_f$snp_pos > lead_snp$snp_pos - window & data_f$snp_pos < lead_snp$snp_pos + window), ]
    # message(paste("Added:", lead_snp$snp_chr, lead_snp$snp_pos))
  }
  return(res)
}


# Declare function definitions
#' Convert Correlation Coefficient to Z-Score
#'
#' This function converts a correlation coefficient `r` to a Z-score based on
#' the given degrees of freedom `df`.
#' The Z-score is useful for comparing analyses on specific significance thresholds.
#'
#' @param r A numeric vector of correlation coefficients.
#' @param df A numeric vector of degrees of freedom corresponding to the correlation coefficients.
#' @return A numeric vector of Z-scores. If the length of `r` and `df` do not match, the function returns 0.
#' @details
#' The function calculates the t-value based on the correlation coefficient and degrees of freedom.
#' It then finds the minimum log p-value from the two-tailed t-distribution and converts it to
#' a Z-score. The sign of the Z-score is adjusted based on the sign of the original correlation coefficient.
#'
#' @examples
#' # Convert a correlation coefficient of 0.5 with 10 degrees of freedom to a Z-score
#' corToZ(0.5, 10)
#'
#' # Convert multiple correlation coefficients to Z-scores
#' corToZ(c(0.3, -0.4), c(15, 20))
#'
#' @export
corToZ <- function(r, df){
  if(length(r) != length(df)){
    return(0L)  #Error condition
  }
  t <- sqrt(df) * abs(r) / sqrt(1 - (r*r))
  logp <- pmin(pt(t, df, log.p = T), pt(t, df, lower.tail=FALSE, log.p = T))#normaly need *2 since this is two-side. However only used to covert to z-score
  z <- qnorm(logp, log.p = T) #no /2 because because logp was not multipled by 2
  z <- ifelse(r > 0, z * -1, z)
  return(z)
}

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


calculate_uncorrected_p_value <- function(cis_r2, df, p_target=5e-8) {
  # Calculate Target z score and r2
  target_z_score <- qnorm(p_target / 2, lower.tail = F)
  target_r2 <- zToCor(rep(target_z_score, length(df)), df)^2

  # Calculate the uncorrected r2
  uncorrected_r2 <- target_r2 / (1 + cis_r2)

  # Calculate the uncorrected z score and p-value
  uncorrected_z_score <- corToZ(sqrt(uncorrected_r2), df)
  uncorrected_p <- 2 * pnorm(q = uncorrected_z_score, lower.tail = FALSE)
  return(list("target_r2" = target_r2, "uncorrected_r2" = uncorrected_r2, "uncorrected_p" = uncorrected_p))
}

correct_for_primary_effect <- function(uncorrected_z_score, primary_r2, df) {
  # Calculate Target z score and r2
  uncorrected_r2 <- zToCor(uncorrected_z_score, df)^2

  # Calculate the uncorrected r2
  corrected_r2 <- uncorrected_r2 * (1 + primary_r2)

  # Calculate the uncorrected z score and p-value
  corrected_z_score <- corToZ(sqrt(corrected_r2), df)
  corrected_p <- 2 * pnorm(q = corrected_z_score, lower.tail = FALSE)
  return(corrected_p)
}

identify_cis_explained_variance <- function(input = "lead_variants.txt") {

  # Process input
  lead_effects <- fread(input)

  # pearson R
  lead_cis_effects <- lead_effects %>%
    filter(type == "cis") %>%
    mutate(
      z_score = beta/se,
      pearson_r = zToCor(z_score, N-1),
      cis_r2 = pearson_r^2) %>%
    rowwise() %>%
    mutate(uncorrected = list(calculate_uncorrected_p_value(cis_r2, df=43301-1, p_target=5e-8))) %>%
    unnest_wider(uncorrected) %>%
    mutate(uncorrected_mlog10p = -log10(uncorrected_p))

  fwrite(lead_cis_effects %>% select(phenotype, cis_r2), "cis_explained_variance_per_gene_20241127.txt", sep="\t", col.names=T, row.names=F, quote=F)

  ggplot(lead_cis_effects, aes(x=uncorrected_mlog10p)) +
    geom_histogram()

}


extract_gene <- function(input_path, gene, cis_explained_variance, variant_reference, gene_reference) {
  uncorrected_p <- calculate_uncorrected_p_value(cis_explained_variance, 43301-1)$uncorrected_p

  trans_boundary <- 5e6

  ds <- open_dataset(sprintf("%s/phenotype=%s", args$input_path, gene))

  extract <- ds %>% collect() %>%
    mutate(z_score = beta/standard_error, p = 2*pnorm(abs(z_score), lower.tail=F)) %>%
    filter(p < uncorrected_p) %>%
    inner_join(variant_reference)

  res <- IdentifyLeadSNPs(
    extract,
    snp_id_col="variant_index",
    snp_chr_col = "chromosome",
    snp_pos_col = "bp",
    beta_col = "beta",
    se_col = "standard_error",
    Pthresh=uncorrected_p)

  results <- extract %>% filter(variant_index %in% res$SNP) %>%
    inner_join(gene_reference) %>%
    mutate(
      cis_like = chromosome == seqname & bp > start - trans_boundary & bp < end + trans_boundary) %>%
    filter(!cis_like) %>%
    mutate(
      corrected_p = correct_for_primary_effect(z_score, cis_explained_variance, sample_size - 1)) %>%
    filter(corrected_p < 5e-8)

  return(results)

}

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }

  args <- parser$parse_args(argv)

  gene_list <- fread(args$genes)$V1

  explained_variance <- fread(args$cis_explained_variance, data.table = F)
  variant_reference <- read_parquet(args$variant_reference)
  gene_reference <- read_parquet(args$gene_reference)

  bound <- bind_rows(mapply(function(gene) {
    cis_explained_variance <- explained_variance[explained_variance$phenotype==gene, "cis_explained_variance"]
    results <- extract_gene(
      args$input, gene, cis_explained_variance,
      variant_reference = variant_reference,
      gene_reference = gene_reference)
  }, gene_list, SIMPLIFY=F))

  fwrite(bound, sprintf("%s.trans_eQTL_list_after_correction.txt.gz", args$output_prefix), sep="\t", col.names=F, row.names=T, quote=F)
  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}