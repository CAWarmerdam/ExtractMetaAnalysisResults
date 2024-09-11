#!/usr/bin/env Rscript


# Load libraries
library(arrow)
library(tidyverse)
library(data.table)
library(ggrastr)

# Declare constants

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
  snp_ref <- read_parquet("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/processed_data/variants/1000G-30x_index.parquet")
  snp_mapping <- snp_ref %>% select(variant_index, variant)

  input_path <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/fine_mapping/validate_in_sample_ld/calculate_ld_emp/output/GTEx_v9_EUR/corr/chr21.locus.tsv"
  dosages <- as.matrix(fread(input_path), rownames="variant")

  AFs <- rowMeans(dosages) / 2
  af_df <- data.frame(
    variant = rownames(dosages),
    AF = AFs,
    mafThreshold = between(AFs, 0.05, 0.95)
  )

  ld_emp <- fread("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/fine_mapping/validate_in_sample_ld/calculate_ld_emp/scripts/empirical.ld")
  ld_emp$V1 <- NULL
  ld_emp_2 <- ld_emp %>% mutate(V1 = colnames(ld_emp)) %>% slice_head(n=1000) %>% pivot_longer(cols = -V1) %>%
    rename(variant1 = "V1", variant2 = "name") %>% filter(variant2 %in% variant1)

  ld_ref <- fread("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/fine_mapping/validate_in_sample_ld/calculate_ld_ref/output/1000G30x/corr/chr21.asis.tsv.gz")
  ld_ref$V1 <- NULL
  ld_ref_2 <- ld_ref %>% mutate(V1 = colnames(ld_emp)) %>% slice_head(n=1000) %>% pivot_longer(cols = -V1) %>%
    rename(variant1 = "V1", variant2 = "name") %>% filter(variant2 %in% variant1)

  ld_un <- fread("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/fine_mapping/validate_in_sample_ld/calculate_ld_emp/output/GTEx_v9_EUR/corr/chr21.asis.tsv.gz")
  ld_un$V1 <- NULL
  ld_un[lower.tri(ld_un)] <- NA
  ld_un_2 <- ld_un %>% mutate(V1 = colnames(ld_un)) %>% slice_head(n=1000) %>% pivot_longer(cols = -V1) %>%
    rename(variant1 = "V1", variant2 = "name") %>% filter(variant2 %in% variant1)

  ld_joined <- ld_ref_2 %>% inner_join(ld_un_2, by =c("variant1", "variant2"), suffix=c("_1000g30x", "_uncorrected"))

  pdf(sprintf("ld_comparison_ref_vs_uncorrected.pdf"))
  p <- ggplot(ld_joined, aes(value_uncorrected^2, value_1000g30x^2)) + rasterize(geom_point(shape=16, alpha=0.2), dpi=300) +
    annotate("text", x=0.9, y=0.1,
             label=sprintf("italic(R) ^ 2 == %.2f",
                           cor(ld_joined$value_uncorrected^2, ld_joined$value_1000g30x^2, use="complete.obs")^2),
             color="#D55E00",
             parse=T)
  plot(p)
  dev.off()

  ld_store_z <- fread("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/fine_mapping/validate_in_sample_ld/calculate_ld_emp/output/GTEx_v9_EUR/chr21.locus.z")
  ld_store <- fread("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/fine_mapping/validate_in_sample_ld/calculate_ld_emp/output/GTEx_v9_EUR/ldstore_analysis/chr21.locus.ld")

  colnames(ld_store) <- ld_store_z$rsid
  ld_store[lower.tri(ld_store)] <- NA
  ld_store_2 <- ld_store %>% mutate(variant1 = ld_store_z$rsid) %>%
    filter(variant1 %in% ld_un_2$variant2) %>%
    select(variant1, any_of(ld_un_2$variant2)) %>%
    pivot_longer(cols = -variant1, values_drop_na=T, names_to='variant2')

  ld_joined <- ld_un_2 %>% inner_join(ld_store_2, by =c("variant1", "variant2"), suffix=c("_uncorrected", "_ldstore"))

  pdf(sprintf("ld_comparison_ldstore_vs_uncorrected.pdf"))
  p <- ggplot(ld_joined, aes(value_uncorrected, value_ldstore)) + rasterize(geom_point(shape=16, alpha=0.2), dpi=300) +
    annotate("text", x=0.9, y=0.1,
             label=sprintf("italic(R) ^ 2 == %.2f",
                           cor(ld_joined$value_uncorrected, ld_joined$value_ldstore, use="complete.obs")^2),
             color="#D55E00",
             parse=T)
  plot(p)
  dev.off()

  ld_files <- Sys.glob("ld_noBatch100ExpPCs*.tsv.gz")

  for (ld_file in ld_files) {

    # Search for the pattern in the filename
    number <- sub("ld_N(\\d+)\\.tsv\\.gz", "\\1", ld_file)
    print(number)

    ld_perm <- fread(ld_file, header=T)
    ld_perm$V1 <- NULL
    ld_perm[lower.tri(ld_perm)] <- NA
    ld_perm_2 <- ld_perm %>% mutate(rowname = as.integer(colnames(ld_perm))) %>%
      inner_join(snp_mapping, by=c("rowname" = "variant_index")) %>%
      rename(variant1 = "variant") %>%
      filter(variant1 %in% ld_emp_2$variant1 | variant1 %in% ld_emp_2$variant2) %>%
      pivot_longer(cols = -variant1, values_drop_na=T) %>%
      mutate(name = as.integer(name)) %>%
      inner_join(snp_mapping, by=c("name" = "variant_index")) %>%
      rename(variant2 = "variant") %>%
      filter(variant2 %in% variant1)

      print(nrow(ld_perm_2))

    ld_joined <- ld_un_2 %>%
      inner_join(ld_perm_2, by =c("variant1", "variant2"), suffix=c("_dosages", "_permuted")) %>%
      inner_join(af_df, by = c("variant1" = "variant")) %>%
      inner_join(af_df, by = c("variant2" = "variant"), suffix=c("", "2")) %>%
      mutate(mafStatus = case_when(mafThreshold & mafThreshold2 ~ "both", mafThreshold | mafThreshold2 ~ "one", TRUE ~ "none"))

    print(nrow(ld_joined))

    ld_joined_summarised <- ld_joined %>% group_by(mafStatus) %>%
      summarise(r2 = cor(value_dosages, value_permuted, use="complete.obs")^2)

    pdf(sprintf("ld_comparison_uncorrected_N%s.pdf", number), width=9, height=4)
    p <- ggplot(ld_joined, aes(value_dosages, value_permuted)) + rasterize(geom_point(shape=16, alpha=0.1, size=0.1), dpi=300) +
      geom_text(data = ld_joined_summarised,
               aes(x=0.8, y=-0.8, label = sprintf("italic(R) ^ 2 == %.3f", r2)),
               color="#D55E00",
               parse=T) +
      coord_fixed() +
      facet_wrap(vars(mafStatus))
    plot(p)
    dev.off()
  }

  # Perform method
  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}