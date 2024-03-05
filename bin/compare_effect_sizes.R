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
library(ggrastr)

# Declare constants
loadfonts()
old <- theme_set(theme_classic(base_size = 10))
theme_update(
  line = element_line(
    colour = "black", size = (0.5 / (ggplot2::.pt * 72.27/96)),
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
  base_path <- "/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/output/empirical_4GenPC20ExpPC_2024-02-27/eqtls"

  master_table <- fread("/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/MetaAnalysis/per_cohort_qc_analysis_empirical_4GenPCNoExpPC_2024-01-29/input/mastertable_empirical_2024-02-22_fixedSampleNames_fixedGeneIds_fixedBiosVersions.txt")
  ds_meta <- open_dataset(file.path(base_path, "meta"))
  ds_cohort <- open_dataset(file.path(base_path, "cohort"))

  z_score_threshold <- qnorm(5e-8/20000/2, lower.tail=FALSE)

  meta_effects <- ds_meta %>%
      mutate(meta_z_score = beta / standard_error) %>%
      select(phenotype, variant, beta, standard_error, meta_z_score) %>%
      filte#!/usr/bin/env Rscript

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
           library(ggrastr)

           # Declare constants
           loadfonts()
           old <- theme_set(theme_classic(base_size = 10))
           theme_update(
             line = element_line(
               colour = "black", size = (0.5 / (ggplot2::.pt * 72.27/96)),
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
             base_path <- "/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/output/empirical_4GenPC20ExpPC_2024-02-27/eqtls"

             master_table <- fread("/gpfs/space/GI/eQTLGen/freeze3/eqtl_mapping/MetaAnalysis/per_cohort_qc_analysis_empirical_4GenPCNoExpPC_2024-01-29/input/mastertable_empirical_2024-02-22_fixedSampleNames_fixedGeneIds_fixedBiosVersions.txt")
             ds_meta <- open_dataset(file.path(base_path, "meta"))
             ds_cohort <- open_dataset(file.path(base_path, "cohort"))

             z_score_threshold <- qnorm(5e-8/20000/2, lower.tail=FALSE)

             meta_effects <- ds_meta %>%
                 mutate(meta_z_score = beta / standard_error) %>%
                 select(phenotype, variant, beta, standard_error, meta_z_score) %>%
                 filter(abs(meta_z_score) > z_score_threshold) %>% collect()

             meta_effects_merged <- ds_cohort %>%
                 select(phenotype, variant, beta, standard_error, cohort) %>%
                 right_join(meta_effects, by=c("phenotype", "variant"), suffix = c("_cohort", "_meta")) %>%
                 mutate(z_score_cohort = beta_cohort / standard_error_cohort) %>%
                 collect() %>%
                 mutate(doe_match = case_when(sign(z_score_cohort) == sign(meta_z_score) ~ "True", sign(z_score_cohort) != sign(meta_z_score) ~ "False", TRUE ~ NA_character_))

             cohort_labels <- meta_effects_merged %>% group_by(cohort) %>%
               summarise(
                   number_of_effects = sum(!is.na(z_score_cohort), na.rm=T),
                   number_of_effects_percent = number_of_effects / n() * 100,
                   concordant = sum(doe_match == "True", na.rm=T),
                   concordant_percent = concordant / number_of_effects * 100,
                   label=sprintf(
                       "nr. of assocs: %s (%.0f%%)\nconcordant assocs: %.0f (%.0f%%)",
                       number_of_effects, number_of_effects_percent,
                       concordant, concordant_percent),
                   x = Inf,
                   y = -Inf,
                   hjust = 1,
                   vjust = 0)

             pdf(sprintf("%s/comparison_z_scores_meta.pdf", "empirical_4GenPC20ExpPC_2024-02-27"), width=9, height=11, useDingbats = F)
             par(xpd = NA)

             ggplot(meta_effects_merged,
                         aes(x=meta_z_score, y=z_score_cohort, colour=doe_match)) +
               geom_point(shape=16, size=0.1, alpha=0.1) +
               scale_colour_manual(values = c("#D55E00", "dimgray"), name = "Matching direction of effect") +
               geom_label(data=cohort_labels,
                   aes(label = label, x=x, y=y),
                   inherit.aes=F,
                   label.padding=unit(0.5, "lines"),
                   label.r=unit(0, "lines"), label.size=0, size=1.8, size.unit='pt',
                   vjust = "inward", hjust = "inward", colour="black", fill=NA) +
               xlab("Z-score (meta-analyzed)") +
               ylab("Z-score (cohort)") +
               facet_wrap(vars(cohort), scales='free') +
               theme(legend.position="bottom")

             dev.off()

           }

           if (sys.nframe() == 0 && !interactive()) {
             main()
           }

r(abs(meta_z_score) > z_score_threshold) %>% collect()

  meta_effects_merged <- ds_cohort %>%
      select(phenotype, variant, beta, standard_error, cohort) %>%
      right_join(meta_effects, by=c("phenotype", "variant"), suffix = c("_cohort", "_meta")) %>%
      mutate(z_score_cohort = beta_cohort / standard_error_cohort) %>%
      collect() %>%
      mutate(doe_match = case_when(sign(z_score_cohort) == sign(meta_z_score) ~ "True", sign(z_score_cohort) != sign(meta_z_score) ~ "False", TRUE ~ NA_character_))

  cohort_labels <- meta_effects_merged %>% group_by(cohort) %>%
    summarise(
        number_of_effects = sum(!is.na(z_score_cohort), na.rm=T),
        number_of_effects_percent = number_of_effects / n() * 100,
        concordant = sum(doe_match == "True", na.rm=T),
        concordant_percent = concordant / number_of_effects * 100,
        label=sprintf(
            "nr. of assocs: %s (%.0f%%)\nconcordant assocs: %.0f (%.0f%%)",
            number_of_effects, number_of_effects_percent,
            concordant, concordant_percent),
        x = Inf,
        y = -Inf,
        hjust = 1,
        vjust = 0)

  pdf(sprintf("%s/comparison_z_scores_meta.pdf", "empirical_4GenPC20ExpPC_2024-02-27"), width=9, height=11, useDingbats = F)
  par(xpd = NA)

  ggplot(meta_effects_merged,
              aes(x=meta_z_score, y=z_score_cohort, colour=doe_match)) +
    geom_point(shape=16, size=0.1, alpha=0.1) +
    scale_colour_manual(values = c("#D55E00", "dimgray"), name = "Matching direction of effect") +
    geom_label(data=cohort_labels,
        aes(label = label, x=x, y=y),
        inherit.aes=F,
        label.padding=unit(0.5, "lines"),
        label.r=unit(0, "lines"), label.size=0, size=1.8, size.unit='pt',
        vjust = "inward", hjust = "inward", colour="black", fill=NA) +
    xlab("Z-score (meta-analyzed)") +
    ylab("Z-score (cohort)") +
    facet_wrap(vars(cohort), scales='free') +
    theme(legend.position="bottom")

  dev.off()

}

if (sys.nframe() == 0 && !interactive()) {
  main()
}

