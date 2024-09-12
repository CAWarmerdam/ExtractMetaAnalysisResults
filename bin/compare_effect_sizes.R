#!/usr/bin/env Rscript

## ----
## Author:  C.A. (Robert) Warmerdam
## Email:   c.a.warmerdam@umcg.nl
##
## Copyright (c) C.A. Warmerdam, 2021
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
library(data.table)
library(tidyverse)
library(arrow)
library(extrafont)
library(rtracklayer)

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
  title = element_text(colour = "#595A5C", face = "bold"),
)


# Main

#' Execute main
#'
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv=NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  base_path <- "/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/output/empirical_qc_per_cohort_4GenPC20ExpPC_2024-06-14/eqtls"

  master_table <- fread("/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/MetaAnalysis/per_cohort_qc_analysis_empirical_4GenPCNoExpPC_2024-01-29/input/mastertable_empirical_2024-02-22_fixedSampleNames_fixedGeneIds_fixedBiosVersions.txt")

  gtf <- "/gpfs/space/GI/eQTLGen/freeze1/InputFilesForPaper/2023-01-28_MetaAnalysis/data/Homo_sapiens.GRCh38.106.gtf.gz"
  gtf <- readGFF(gtf)
  gtf2 <- gtf[gtf$type == "gene", ]
  gtf2 <- unique(gtf2[, c(1, 4, 5, 9, 11)])

  snp_ref <- open_dataset("/gpfs/space/GI/eQTLGen/processed_data/variants/1000G-30x.index.gz") %>% collect()
  ds_meta <- open_dataset(file.path(base_path, "meta"))
  ds_cohort <- open_dataset(file.path(base_path, "cohort"))
  ds_qgs <- open_dataset(file.path("/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/output/empirical_qc_per_cohort_4GenPC20ExpPC_2024-06-20/eqtls", "cohort"))

  z_score_threshold <- qnorm(5e-8/20000/2, lower.tail=FALSE)

  meta_effects <- ds_meta %>%
      mutate(meta_z_score = beta / standard_error) %>%
      select(phenotype, variant_index, beta, standard_error, meta_z_score, i_squared) %>%
      filter(abs(meta_z_score) > z_score_threshold) %>% collect()

  meta_effects_merged <- ds_cohort %>%
    select(phenotype, variant_index, beta, standard_error, cohort) %>%
    mutate(variant_index = as.integer(variant_index)) %>%
    right_join(meta_effects, by=c("phenotype", "variant_index"), suffix = c("_cohort", "_meta")) %>%
    mutate(z_score_cohort = beta_cohort / standard_error_cohort) %>%
    collect() %>%
    mutate(doe_match = case_when(sign(z_score_cohort) == sign(meta_z_score) ~ "True", sign(z_score_cohort) != sign(meta_z_score) ~ "False", TRUE ~ NA_character_)) %>%
    inner_join(snp_ref)

  meta_effects_merged_2 <- df_qgs %>% filter(cohort == "QGP") %>% mutate(cohort = "QGP3") %>% bind_rows(meta_effects_merged) %>%
    mutate(include_in_meta_analysis = !(cohort %in% c("QGP"))) %>%
    group_by(phenotype, variant_index) %>%
    mutate(weight = 1/(standard_error_cohort^2), beta_weighted = weight * beta_cohort) %>%
    mutate(beta_meta_2 = sum(if_else(include_in_meta_analysis, beta_weighted, 0),na.rm=T) / sum(if_else(include_in_meta_analysis, weight, 0),na.rm=T),
           se_meta_2 = sqrt(1/sum(if_else(include_in_meta_analysis, weight, 0),na.rm=T)),
           meta_z_score = beta_meta_2 / se_meta_2) %>%
    inner_join(gtf2, by = c("phenotype" = "gene_id")) %>%
    mutate(type = case_when(chromosome == seqid & between(bp, start - 1e6, end + 1e6) ~ "cis", TRUE ~ "trans")) %>%
    filter(abs(meta_z_score) > z_score_threshold)

  meta_effects_merged_2 <- meta_effects_merged %>%
    mutate(include_in_meta_analysis = !(cohort %in% c("QGP", "QGP2"))) %>%
    group_by(phenotype, variant_index) %>%
    mutate(weight = 1/(standard_error_cohort^2), beta_weighted = weight * beta_cohort) %>%
    mutate(beta_meta_2 = sum(if_else(include_in_meta_analysis, beta_weighted, 0),na.rm=T) / sum(if_else(include_in_meta_analysis, weight, 0),na.rm=T),
              se_meta_2 = sqrt(1/sum(if_else(include_in_meta_analysis, weight, 0),na.rm=T)),
           meta_z_score = beta_meta_2 / se_meta_2)

  cohort_labels <- meta_effects_merged_2 %>% group_by(cohort, type) %>%
    summarise(
        number_of_effects = sum(!is.na(z_score_cohort), na.rm=T),
        number_of_effects_percent = number_of_effects / n() * 100,
        concordant = sum(doe_match == "True", na.rm=T),
        concordant_percent = concordant / number_of_effects * 100,
        slope = lm(beta_cohort ~ beta_meta)$coefficients["beta_meta"],
        label=sprintf(
            "nr. of assocs: %s (%.0f%%)\nconcordant assocs: %.0f (%.0f%%)\nslope of betas: %.2f",
            number_of_effects, number_of_effects_percent,
            concordant, concordant_percent, slope),
        x = Inf,
        y = -Inf,
        hjust = 1,
        vjust = 0)

  pdf(sprintf("%s/comparison_z_scores_meta.pdf", "empirical_4GenPC20ExpPC_2024-02-27"), width=9, height=11, useDingbats = F)
  par(xpd = NA)

  png(sprintf("%s/comparison_z_scores_meta_cis_5.png", "."), units = 'in', res = 300, width=9, height=11, type='cairo')

  ggplot(meta_effects_merged_2 %>% filter(type == "cis"),
              aes(x=meta_z_score, y=z_score_cohort, colour=doe_match)) +
    geom_point(shape=16, size=0.1, alpha=0.1) +
    scale_colour_manual(values = c("#D55E00", "dimgray"), name = "Matching direction of effect") +
    geom_label(data=cohort_labels %>% filter(type == "cis"),
        aes(label = label, x=x, y=y),
        inherit.aes=F,
        label.padding=unit(0.5, "lines"),
        label.r=unit(0, "lines"), label.size=0, size=2, size.unit='pt',
        vjust = "inward", hjust = "inward", colour="black", fill=NA) +
    ggtitle("Comparison of cis-effects across cohorts") +
    xlab("Z-score (meta-analyzed)") +
    ylab("Z-score (cohort)") +
    facet_wrap(vars(cohort), scales='free') +
    theme(legend.position="bottom")

  dev.off()

  png(sprintf("%s/comparison_z_scores_meta_trans_5.png", "."), units = 'in', res = 300, width=9, height=11, type='cairo')

  ggplot(meta_effects_merged_2 %>% filter(type == "trans"),
         aes(x=meta_z_score, y=z_score_cohort, colour=doe_match)) +
    geom_point(shape=16, size=0.1, alpha=0.1) +
    scale_colour_manual(values = c("#D55E00", "dimgray"), name = "Matching direction of effect") +
    geom_label(data=cohort_labels %>% filter(type == "trans"),
               aes(label = label, x=x, y=y),
               inherit.aes=F,
               label.padding=unit(0.5, "lines"),
               label.r=unit(0, "lines"), label.size=0, size=2, size.unit='pt',
               vjust = "inward", hjust = "inward", colour="black", fill=NA) +
    ggtitle("Comparison of trans-effects across cohorts") +
    xlab("Z-score (meta-analyzed)") +
    ylab("Z-score (cohort)") +
    facet_wrap(vars(cohort), scales='free') +
    theme(legend.position="bottom")

  dev.off()


}

if (sys.nframe() == 0 && !interactive()) {
  main()
}

