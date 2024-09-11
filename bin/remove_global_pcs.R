#!/usr/bin/env Rscript


# Load libraries
library(argparse)
library(tidyverse)
library(data.table)
library(dtplyr)
library(arrow)
library(Rfast)
library(rtracklayer)

# Declare constants
BASE_DIR <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2"
SNP_REF <- file.path(BASE_DIR, "processed_data/variants/1000G-30x_index.parquet")
GENE_REF <- file.path(BASE_DIR, "public_data/Homo_sapiens.GRCh38.106.gtf.gz")

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
  logp <- pt(abs(t), df, lower.tail = FALSE, log.p = T) #normaly need *2 since this is two-side. However only used to covert to z-score
  #logp <- pmin(pt(t, df, log.p = T), pt(t, df, lower.tail=FALSE, log.p = T)) #normaly need *2 since this is two-side. However only used to covert to z-score
  z <- qnorm(logp, log.p = T) #no /2 because because logp was not multipled by 2
  z <- z * -sign(r)
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

#' Correct Post-Hoc eQTL Z-Scores Using Principal Components
#'
#' This function applies a post-hoc correction to eQTL Z-scores by accounting for
#' the effects of principal components (PCs) derived from a gene-gene correlation matrix.
#' The function corrects the Pearson correlations by removing the contribution of the PCs,
#' recalculates the Z-scores, and returns a data.table with corrected and partially corrected
#' Z-scores, along with other relevant information.
#'
#' @param pearson_correlations A matrix of eQTL Z-scores converted to Pearson correlations.
#'                             The matrix should have genes as rows and variants as columns.
#' @param eigenvectors A matrix of eigenvectors of the gene-gene correlation matrix.
#'                     Each column corresponds to a principal component (PC).
#' @param sample_size_matrix A matrix indicating the sample sizes for each eQTL,
#'                           with genes as rows and variants as columns.
#' @param cis_trans_matrix A binary matrix indicating whether each eQTL is cis (TRUE) or trans (FALSE).
#'                         The matrix should have genes as rows and variants as columns.
#' @param number_of_pcs An integer specifying the number of principal components to include in the correction.
#' @param component_variant_estimates A matrix representing the effects (slopes) of the principal components on each SNP.
#'               The matrix should have PCs as rows and SNPs as columns.
#'
#' @return A data.table containing the following columns:
#'   - `row_names`: The row names of the original matrices (genes).
#'   - `col_names`: The column names of the original matrices (variants).
#'   - `sample_size`: The sample size for each eQTL.
#'   - `cis`: A boolean indicating whether the eQTL is cis (TRUE) or trans (FALSE).
#'   - `z_scores_corrected`: The fully corrected Z-scores.
#'   - `z_scores_corrected_partial`: The partially corrected Z-scores.
#'
#' @export
correct_post_hoc <- function(eqtls_as_corr_matrix, eigenvectors, sample_size_matrix, cis_trans_matrix,
                             number_of_pcs, component_variant_estimates, unbiased_snps, test_snps=NULL) {

  # Calculate the contribution of the PCs to the eQTLs by multiplying the eigenvectors
  # with the slopes for the specified number of PCs.
  eqtls_predicted <- t(eigenvectors[colnames(eqtls_as_corr_matrix), 1:number_of_pcs] %*% component_variant_estimates[1:number_of_pcs,])

  # Partially correct the Pearson correlations by subtracting the PC contributions.
  if (!is.null(test_snps)) {
    eqtls_corr_corrected_partial <- eqtls_as_corr_matrix[test_snps, ] - eqtls_predicted[test_snps, ]
  } else {
    eqtls_corr_corrected_partial <- eqtls_as_corr_matrix - eqtls_predicted
  }

  # Calculate the residual variance explained after accounting for the PCs.
  residual_variance_explained <- 1 - corpairs(
    eqtls_as_corr_matrix[unbiased_snps, ],
    eqtls_predicted[unbiased_snps, ])^2
  names(residual_variance_explained) <- colnames(eqtls_predicted)

  # eqtls_corr_corrected_bad <- sign(eqtls_corr_corrected_partial) * sqrt(eqtls_corr_corrected_partial^2 / residual_variance_explained)
  # Fully correct the Pearson correlations by adjusting for residual variance.
  eqtls_corr_corrected <- sign(eqtls_corr_corrected_partial) * t(sqrt(t(eqtls_corr_corrected_partial^2) / residual_variance_explained))
  # # Fully correct the Pearson correlations by adjusting for residual variance snp wise

  # Clamp the corrected correlations to the range [-1, 1] to ensure valid correlation values.
  eqtls_corr_corrected <- pmin(pmax(eqtls_corr_corrected, -1), 1)

  # Extract row and column names for the final data.table.
  output_rows <- rownames(eqtls_corr_corrected)
  output_cols <- colnames(eqtls_corr_corrected)

  # Convert the fully corrected correlations back to Z-scores.
  eqtl_corrected_z_scores <- corToZ(
    eqtls_corr_corrected[output_rows, output_cols],
    sample_size_matrix[output_rows, output_cols] - 1)

  rm(eqtls_corr_corrected)

  # Convert the partially corrected correlations back to Z-scores.
  eqtl_partially_corrected_z_scores <- corToZ(
    eqtls_corr_corrected_partial[output_rows, output_cols],
    sample_size_matrix[output_rows, output_cols] - 1)

  rm(eqtls_corr_corrected_partial)

  # Convert the uncorrected correlations back to Z-scores.
  eqtls_uncorrected_z_scores <- corToZ(
    eqtls_as_corr_matrix[output_rows, output_cols] ,
    sample_size_matrix[output_rows, output_cols] - 1)

  rm(eqtls_as_corr_matrix)

  # Convert the uncorrected correlations back to Z-scores.
  predicted_eqtls_as_z_scores <- corToZ(
    eqtls_predicted[output_rows, output_cols],
    sample_size_matrix[output_rows, output_cols] - 1)

  rm(eqtls_predicted)

  # Combine the results into a data.table.
  result_dt <- cbind(
    as.data.table(
      expand.grid(variant_index = output_rows, phenotype = output_cols)
    ),
    sample_size = as.vector(sample_size_matrix[output_rows, output_cols]),  # Include the sample sizes.
    cis = as.vector(cis_trans_matrix[output_rows, output_cols]),  # Indicate cis/trans status.
    z_scores_uncorrected = as.vector(eqtls_uncorrected_z_scores[output_rows, output_cols]), # original
    z_scores_predicted = as.vector(predicted_eqtls_as_z_scores[output_rows, output_cols]), # he ho
    z_scores_corrected = as.vector(eqtl_corrected_z_scores[output_rows, output_cols]),  # Fully corrected Z-scores.
    z_scores_corrected_partial = as.vector(eqtl_partially_corrected_z_scores[output_rows, output_cols])  # Partially corrected Z-scores.
  )

  result_dt$residual_variance_explained <- residual_variance_explained[result_dt$phenotype]

  return(result_dt)
}

#' Maximum Independent Set of Uncorrelated Genes
#'
#' This function identifies the maximum set of uncorrelated genes (independent vertices)
#' from a boolean adjacency matrix. It uses an algorithm that selects the vertex with the
#' least amount of edges and removes all vertices that are connected to it.
#'
#' @param adjacency_matrix A boolean adjacency matrix representing gene correlations (1 if correlated, 0 if not).
#' @param nodes A character vector containing the names of the genes corresponding to the rows/columns of the matrix.
#' @return A character vector of genes that are part of the maximum independent set.
#' @examples
#' adj_matrix <- matrix(c(1, 1, 0, 1, 1, 1, 0, 0, 1), nrow = 3, byrow = TRUE)
#' gene_names <- c("Gene1", "Gene2", "Gene3")
#' maximum_independent_set(adj_matrix, gene_names)
maximum_independent_set <- function(adjacency_matrix, nodes) {
  # Initialize the set of independent vertices to be empty
  independent_vertices <- character(0)
  active_nodes <- nodes == nodes
  names(active_nodes) <- nodes

  colnames(adjacency_matrix) <- nodes
  rownames(adjacency_matrix) <- nodes

  diag(adjacency_matrix) <- T

  # First collect every gene that is totally independent.
  # Removing these will not have an effect on other genes.
  unconnected_vertices <- rowSums(adjacency_matrix) == 1
  independent_vertices <- c(independent_vertices, nodes[unconnected_vertices])

  with_edges <- colSums(adjacency_matrix[unconnected_vertices, active_nodes] == T) > 0
  active_nodes[with_edges] <- F

  print(length(unconnected_vertices))
  print(sum(active_nodes))

  # Loop while the matrix is not empty
  while (any(active_nodes)) {
    # Find the vertex with the least amount of edges
    least_connected <- names(base::which.min(colSums(adjacency_matrix[active_nodes, active_nodes, drop=F])))
    cat(least_connected, "\n")

    # Add the least connected vertex to the independent set
    independent_vertices <- c(independent_vertices, least_connected)

    # Remove all vertices that have an edge with the least connected vertex
    with_edges <- adjacency_matrix[least_connected,][active_nodes, drop=F] == T

    active_nodes[names(with_edges)[with_edges]] <- F

    cat(sum(active_nodes), "\n")
  }

  return(independent_vertices)
}


calculate_component_variant_estimates_lco <- function(
  eqtls_as_corr_matrix, leave_chromosome_out, eigenvectors_with_intercept, uncorrelated_genes, eigen_subset_func) {

  component_variant_estimates <- matrix(NA, ncol = nrow(eqtls_as_corr_matrix), nrow = ncol(eigenvectors),
                                        dimnames = list(colnames(eigenvectors), rownames(eqtls_as_corr_matrix)))

  all_genes <- colnames(eqtls_as_corr_matrix)
  all_genes_uncorrelated <- all_genes[all_genes %in% uncorrelated_genes]

  for (i in seq_len(nrow(leave_chromosome_out))) {
    variants <- as.character(unlist(leave_chromosome_out[i, 'variant_index']))
    leave_genes_out <- unname(unlist(leave_chromosome_out[i, 'gene_id']))

    # Identify genes to exclude
    genes_to_include <- all_genes_uncorrelated[!all_genes_uncorrelated %in% leave_genes_out]

    # Subset matrices
    eigen_subset <- eigenvectors_with_intercept[genes_to_include,]
    eqtl_subset <- eqtls_as_corr_matrix[variants, genes_to_include, drop = FALSE]

    # Calculate component estimates for the SNP
    if (ncol(eigen_subset) > 0) {  # Ensure there are genes to include
      component_variant_estimates[,variants] <- eigen_subset_func(t(eqtl_subset), eigen_subset)
    } else {
      component_variant_estimates[,variants] <- NA  # No genes to include, result is NA
    }

    # Print the progress (only for "corr" estimate)
    print(paste("Processed", i, "out of", nrow(leave_chromosome_out), "chromosomes"))
  }

  # Save the results
  return(component_variant_estimates)
}


calculate_component_variant_estimates <- function(
  eqtls_as_corr_matrix, cis_trans_matrix, eigenvectors_with_intercept, uncorrelated_genes, eigen_subset_func) {
  component_variant_estimates <- matrix(NA, ncol = nrow(eqtls_as_corr_matrix), nrow = ncol(eigenvectors),
                                        dimnames = list(colnames(eigenvectors), rownames(eqtls_as_corr_matrix)))

  for (i in seq_len(nrow(eqtls_as_corr_matrix))) {
    snp <- rownames(eqtls_as_corr_matrix)[i]

    # Identify genes to exclude
    genes_without_cis <- colnames(cis_trans_matrix)[cis_trans_matrix[i, colnames(eqtls_as_corr_matrix)] == 0]
    genes_to_include <- genes_without_cis[genes_without_cis %in% uncorrelated_genes]

    # Subset matrices
    eigen_subset <- eigenvectors_with_intercept[genes_to_include,]
    eqtl_subset <- as.vector(eqtls_as_corr_matrix[i, genes_to_include, drop = FALSE])

    # Calculate component estimates for the SNP
    if (ncol(eigen_subset) > 0) {  # Ensure there are genes to include
      component_variant_estimates[,snp] <- eigen_subset_func(eqtl_subset, eigen_subset)
    } else {
      component_variant_estimates[,snp] <- NA  # No genes to include, result is NA
    }

    # Print the progress (only for "corr" estimate)
    if (i %% 100 == 0) {
      print(paste("Processed", i, "out of", nrow(eqtls_as_corr_matrix), "variants"))
    }
  }

  # Save the results
  return(component_variant_estimates)
}

# Function to read the files based on input prefix
read_from_prefix <- function(prefix) {
  # Define file paths
  mat_n_txt_path <- paste0(prefix, ".n.txt.gz")
  mat_z_txt_path <- paste0(prefix, ".z.txt.gz")
  mat_n_h5_path <- paste0(prefix, ".n.h5")
  mat_z_h5_path <- paste0(prefix, ".z.h5")

  # Initialize matrices
  mat_n <- NULL
  mat_z <- NULL

  # Read n matrix
  if (file.exists(mat_n_h5_path)) {
    mat_n <- h5read(mat_n_h5_path, "sample_size")
  } else if (file.exists(mat_n_txt_path)) {
    mat_n <- fread(mat_n_txt_path)
  }

  # Read z matrix
  if (file.exists(mat_z_h5_path)) {
    mat_z <- h5read(mat_z_h5_path, "z_scores")
  } else if (file.exists(mat_z_txt_path)) {
    mat_z <- fread(mat_z_txt_path)
  }

  return(list(mat_n = mat_n, mat_z = mat_z))
}

# Main

#' Execute main
#' 
#' @param argv A vector of arguments normally supplied via command-line.
main <- function(argv = NULL) {
  if (is.null(argv)) {
    argv <- commandArgs(trailingOnly = T)
  }
  # Process input

  # parse command-line arguments
  parser <- argparse::ArgumentParser(description='Calculate correlations between z-scores and find uncorrelated genes')
  parser$add_argument('--input-prefix', dest='input_prefix', required=FALSE,
                      help='Prefix for input files (e.g., "mat").')
  parser$add_argument('--input-parquet', dest='input_parquet', required=FALSE,
                      help='Path to input Parquet file.')
  parser$add_argument('--variant-ref', dest='variant_reference', required=TRUE)
  parser$add_argument('--gene-ref', dest='gene_reference', required=TRUE)
  parser$add_argument('--output-file', dest='output_file', help='Path to the output file')
  parser$add_argument('--gene-pca-prefix', dest='components', help='Path to gene components',
                      default=None, required=TRUE)
  parser$add_argument('-c', '--n-components', type=float, default=100, help='Number of components to regress out')
  parser$add_argument('-n', '--n-threshold', type=float, default=0, help='Minimum sample size of gene-variant pairs to consider for gene-gene correlations')
  args <- parser$parse_args()

  estimate <- "beta"

  uncorrelated_genes_file_path <- "uncorrelated_genes_r2Max0.05.rds"
  if (file.exists(uncorrelated_genes_file_path)) {
    uncorrelated_genes <- readRDS(uncorrelated_genes_file_path)
  } else {
    corr_matrix <- readRDS(file.path(BASE_DIR, "freeze3/Interpretation/extractions/correction_optimization/permuted1_RNAseq_4GenPCTechCov_2024-06-26/exported_matrix/gene_correlation_matrix.rds"))
    uncorrelated_genes <- maximum_independent_set(abs(corr_matrix) > sqrt(0.05), rownames(corr_matrix))
    saveRDS(uncorrelated_genes, uncorrelated_genes_file_path)
  }

  eigenvectors <- readRDS("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/correction_optimization/permuted1_RNAseq_4GenPCTechCov_2024-06-26/exported_matrix/gene_eigenvectors.rds")
  colnames(eigenvectors) <- paste0("component", c(1:ncol(eigenvectors)))
  eigenvectors_with_intercept <- cbind(intercept=1, eigenvectors)  # Add a column of 1s for the intercept

  #eigen <- readRDS(args$components)

  eqtls <- open_dataset(file.path(BASE_DIR, "freeze3/eqtl_mapping/output/empirical/empirical0_RNAseq_4GenPCTechCov_2024-08-21/eqtls/meta/")) %>% mutate(z_score = beta/standard_error) %>% collect() %>% as.data.table()

  z_score_matrix <- as.matrix(dcast(eqtls, variant_index ~ phenotype, value.var = "z_score"), rownames=1)
  sample_size_matrix <- as.matrix(dcast(eqtls, variant_index ~ phenotype, value.var = "sample_size"), rownames=1)
  rm(eqtls)

  variants_to_load <- rownames(z_score_matrix)
  snp_ref <- read_parquet(SNP_REF) %>% filter(variant_index %in% variants_to_load)
  snp_ref_dt <- snp_ref %>% mutate(snp_start = bp, snp_end = bp, chromosome=as.character(chromosome)) %>% as.data.table()

  gene_ref <- readGFF(GENE_REF)
  gene_ref <- gene_ref[gene_ref$type == "gene", ]
  gene_ref <- unique(gene_ref[, c(1, 4, 5, 9, 11)])
  gene_ref_dt <- gene_ref %>%
    filter(gene_id %in% colnames(z_score_matrix)) %>%
    mutate(chromosome = as.character(seqid),
           trans_window_upper = end + 5e6,
           trans_window_lower = start - 5e6) %>%
    data.table()

  setkey(snp_ref_dt, chromosome)
  setkey(gene_ref_dt, chromosome)

  leave_chromosome_out <- merge(snp_ref_dt, gene_ref_dt, allow.cartesian = T) %>%
    group_by(chromosome) %>%
    summarise(variant_index = list(unique(variant_index)), gene_id = list(unique(gene_id))) %>% collect()

  z_score_matrix <- z_score_matrix[,colnames(z_score_matrix) %in% gene_ref_dt$gene_id]
  sample_size_matrix <- sample_size_matrix[,colnames(sample_size_matrix) %in% gene_ref_dt$gene_id]

  setkey(snp_ref_dt, chromosome, snp_start, snp_end)
  setkey(gene_ref_dt, chromosome, trans_window_lower, trans_window_upper)

  overlap_result <- foverlaps(snp_ref_dt, gene_ref_dt, nomatch = 0)

  overlap_result_mat <- as.matrix(dcast(overlap_result, variant_index ~ gene_id,
                              value.var = "gene_id",
                              fun.aggregate = length), rownames=1)

  cis_trans_matrix <- matrix(FALSE, nrow = nrow(z_score_matrix), ncol = ncol(z_score_matrix),
                             dimnames = list(rownames(z_score_matrix), colnames(z_score_matrix)))

  cis_trans_matrix[rownames(overlap_result_mat), colnames(overlap_result_mat)] <- as.logical(overlap_result_mat)

  eqtls_as_corr_matrix <- zToCor(z_score_matrix, sample_size_matrix - 1)

  if (estimate == "beta") {
    eigen_subset_func <- function(eqtl_subset, eigen_subset) {
      beta_hat <- solve(t(eigen_subset) %*% eigen_subset) %*% t(eigen_subset) %*% eqtl_subset
      return(beta_hat[-1,])  # Remove the first row corresponding to the intercept
    }
    save_path <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/beta_hat_matrix_all_genes.rds"
    save_path <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/beta_hat_matrix.rds"
  } else if (estimate == "corr") {
    eigen_subset_func <- function(eqtl_subset, eigen_subset) {
      correlation_coefficients <- cor(eqtl_subset, eigen_subset[,-1])
      return(correlation_coefficients)
    }
    save_path <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/correlations_matrix_large.rds"
  } else if (estimate == "cov") {
    eigen_subset_func <- function(eqtl_subset, eigen_subset) {
      covariance_estimates <- cov(eqtl_subset, eigen_subset[,-1])
      return(covariance_estimates)
    }
    save_path <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/covariance_matrix.rds"
  }

  if (file.exists(save_path)) {
    component_variant_estimates <- readRDS(save_path)
  } else {
    # component_variant_estimates <- calculate_component_variant_estimates(
    #   eqtls_as_corr_matrix, cis_trans_matrix, eigenvectors_with_intercept, uncorrelated_genes, eigen_subset_func)

    component_variant_estimates <- calculate_component_variant_estimates_lco(
      eqtls_as_corr_matrix, leave_chromosome_out, eigenvectors_with_intercept, uncorrelated_genes, eigen_subset_func)

    saveRDS(component_variant_estimates, save_path)
  }

  # Define a vector of different numbers of PCs to correct for
  pcs_vector <- c(20, 50, 100, 250, 500, 1000)  # Replace with your desired values

  unbiased_snps <- fread("/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/correction_optimization/permuted1_RNAseq_4GenPCTechCov_2024-06-26/exported_matrix/uncorrelated_variants_passing_sample_size_threshold.tsv.gz", col.names=c("variant")) %>%
    inner_join(snp_ref) %>%
    filter(variant_index %in% rownames(eqtls_as_corr_matrix)) %>%
    pull(variant_index)

  base_dir <- "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2"
  pruned_immune_snps <- fread(file.path(base_dir, "correction_optimization/pruned_immune_snps/immune_snps_chrAll.prune.in"),
                       header=F, col.names = c("variant")) %>%
    inner_join(snp_ref) %>%
    filter(variant_index %in% rownames(eqtls_as_corr_matrix)) %>%
    pull(variant_index)


  results <- rbindlist(
    lapply(pcs_vector, function(n_pcs) {
      print(message(n_pcs))
      post_hoc_out <- correct_post_hoc(eqtls_as_corr_matrix,
                                 eigenvectors,
                                 sample_size_matrix,
                                 cis_trans_matrix,
                                 n_pcs,
                                 component_variant_estimates,
                                 unbiased_snps = as.character(unbiased_snps),
                                 test_snps = as.character(pruned_immune_snps))
      post_hoc_out[, corrected_pcs := n_pcs]  # Add a column indicating the number of PCs corrected
      return(post_hoc_out)
    })
  )

  saveRDS(results, "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/eqtls_post_hoc_corrected_comparison_uncorrelated_genes_independent_snps_beta.rds")
  fwrite(results, "/groups/umcg-fg/tmp01/projects/eqtlgen-phase2/analyses/correction_optimization/eqtls_post_hoc_corrected_comparison_uncorrelated_genes_independent_snps_beta.txt.gz",
         row.names=T, col.names=T, quote=F, sep="\t")

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}