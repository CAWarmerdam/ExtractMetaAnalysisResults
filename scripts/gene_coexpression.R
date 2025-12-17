#!/usr/bin/env Rscript


# Load libraries
library(tidyverse)
library(data.table)
library(arrow)

# Declare constants

# Declare function definitions


pca_train <- function(correlation_matrix = NULL,
                      max_nr_comps = NULL,
                      variance_threshold = NULL) {
  # Perform Eigen decomposition on the correlation matrix
  eigen_decomp <- eigen(correlation_matrix)
  eigenvalues <- eigen_decomp$values
  eigenvectors <- eigen_decomp$vectors
  rownames(eigenvectors) <- rownames(correlation_matrix)
  names(eigenvalues) <- rownames(correlation_matrix)

  # Calculate cumulative explained variance
  explained_variance <- eigenvalues / sum(eigenvalues)
  cumulative_explained_variance <- cumsum(explained_variance)

  # Determine the number of components to retain
  if (!is.null(variance_threshold)) {
    # Find the number of components required to meet the threshold
    nr_comps <- which(cumulative_explained_variance >= variance_threshold)[1]
  } else if (!is.null(max_nr_comps)) {
    # Limit by max number of components if specified
    nr_comps <- min(max_nr_comps, ncol(correlation_matrix))
  } else {
    # Default to all components
    nr_comps <- ncol(correlation_matrix)
  }

  # Keep only the top `nr_comps` components
  eigenvectors_reduced <- eigenvectors[, 1:nr_comps]
  eigenvalues_reduced <- eigenvalues[1:nr_comps]

  return(list(
    eigenvalues = eigenvalues_reduced,
    explained_variance = explained_variance,
    cumulative_explained_variance = cumulative_explained_variance,
    eigenvectors = eigenvectors_reduced,
    retained_components = nr_comps
  ))
}

svd_train <- function(correlation_matrix = NULL,
                      max_nr_comps = NULL,
                      variance_threshold = NULL) {

  # Perform Eigen decomposition on the correlation matrix
  svd_out <- svd(correlation_matrix)
  eigenvalues <- svd_out$d
  eigenvectors <- svd_out$u
  rownames(eigenvectors) <- rownames(correlation_matrix)
  names(eigenvalues) <- rownames(correlation_matrix)

  # Calculate cumulative explained variance
  explained_variance <- eigenvalues / sum(eigenvalues)
  cumulative_explained_variance <- cumsum(explained_variance)

  # Determine the number of components to retain
  if (!is.null(variance_threshold)) {
    # Find the number of components required to meet the threshold
    nr_comps <- which(cumulative_explained_variance >= variance_threshold)[1]
  } else if (!is.null(max_nr_comps)) {
    # Limit by max number of components if specified
    nr_comps <- min(max_nr_comps, ncol(correlation_matrix))
  } else {
    # Default to all components
    nr_comps <- ncol(correlation_matrix)
  }

  # Keep only the top `nr_comps` components
  eigenvectors_reduced <- eigenvectors[, 1:nr_comps]
  eigenvalues_reduced <- eigenvalues[1:nr_comps]

  return(list(
    eigenvalues = eigenvalues_reduced,
    explained_variance = explained_variance,
    cumulative_explained_variance = cumulative_explained_variance,
    eigenvectors = eigenvectors_reduced,
    retained_components = nr_comps
  ))
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

# T to cor
tToCor <- function(t, df){
  if(length(t) != length(df)){
    return(0L)  #Error condition
  }
  r <- t/sqrt(df+t^2)
  return(r)
}


write_pca_fit <- function(eigenvectors, eigenvalues, gene_centers, gene_deviations, out_label) {
  genes <- rownames(eigenvectors)

  component_names <- paste0("comp", seq_len(ncol(eigenvectors)))
  colnames(eigenvectors) <- component_names
  eigenvectors_dt <- as.data.table(eigenvectors, keep.rownames="gene")
  names(eigenvalues) <- component_names
  eigenvalues_dt <- as.data.table(eigenvalues, keep.rownames = "component")
  colnames(eigenvalues_dt) <- c("component", "value")

  str(names(eigenvalues_dt))

  gene_centers_dt <- as.data.table(gene_centers[genes], keep.rownames = "gene")
  colnames(gene_centers_dt) <- c("gene", "value")
  gene_deviations_dt <- as.data.table(gene_deviations[genes], keep.rownames = "gene")
  colnames(gene_deviations_dt) <- c("gene", "value")

  fwrite(eigenvectors_dt, sprintf("%s_eigenvectors.txt.gz", out_label), sep="\t", row.names=F)
  fwrite(eigenvalues_dt, sprintf("%s_eigenvalues.txt.gz", out_label), sep="\t", row.names=F)
  fwrite(gene_centers_dt, sprintf("%s_gene_centers.txt.gz", out_label), sep="\t", row.names=F)
  fwrite(gene_deviations_dt, sprintf("%s_gene_deviations.txt.gz", out_label), sep="\t", row.names=F)

  write_feather(eigenvectors_dt, sprintf("%s_eigenvectors.feather", out_label))
  write_feather(eigenvalues_dt, sprintf("%s_eigenvalues.feather", out_label))
  write_feather(gene_centers_dt, sprintf("%s_gene_centers.feather", out_label))
  write_feather(gene_deviations_dt, sprintf("%s_gene_deviations.feather", out_label))
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
  mat_n <- fread("mat.n.txt.gz") %>% as.matrix(rownames=1)
  mat_z <- fread("mat.z.txt.gz") %>% as.matrix(rownames=1)
  all(colnames(mat_n) == colnames(mat_z))
  all(rownames(mat_n) == rownames(mat_z))

  # Perform method

  rho_mat <- tToCor(mat_z, mat_n-1)

  # Center and scale each gene
  rho_mat <- scale(rho_mat)
  rho_mat_attributes <- attributes(rho_mat)

  # Calculate correlations
  coexp_all_genes <- crossprod(rho_mat) / (nrow(rho_mat) - 1)

  saveRDS(rho_mat_attributes$`scaled:center`, "gene_centers_16781_genes.rds")
  saveRDS(rho_mat_attributes$`scaled:scale`, "gene_deviation_16781_genes.rds")

  gene_average_n <- colMeans(mat_n)
  genes_pass_10k <- names(gene_average_n[gene_average_n > (43301 * 0.90)])

  saveRDS(coexp_all_genes, "coexpression_matrix_16781_genes.rds")
  saveRDS(coexp_all_genes[genes_pass_10k, genes_pass_10k], "coexpression_matrix_10497_genes.rds")

  saveRDS(rho_mat_attributes$`scaled:center`[genes_pass_10k], "gene_centers_10497_genes.rds")
  saveRDS(rho_mat_attributes$`scaled:scale`[genes_pass_10k], "gene_deviation_10497_genes.rds")

  rm(mat_n, mat_z)

  # Process output
  coexp_dt <- data.table(gene = rownames(coexp_all_genes), coexp_all_genes)
  coexp_10k_dt <- data.table(gene = genes_pass_10k, coexp_all_genes[genes_pass_10k, genes_pass_10k])

  fwrite(coexp_dt, "coexpression_matrix_16781_genes.txt.gz", sep="\t", row.names=F, col.names=T, quote=F)
  fwrite(coexp_10k_dt, "coexpression_matrix_10497_genes.txt.gz", sep="\t", row.names=F, col.names=T, quote=F)

  rm(list=ls())

  coexp_all_genes <- readRDS( "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/all_empirical_independentVariants_4GenPC20ExpPC_2025-01-22/exported_matrix/coexpression_matrix_16781_genes.rds")
  coexp_10k <- readRDS( "/scratch/hb-functionalgenomics/projects/eqtlgen-phase2/freeze3/Interpretation/extractions/all_empirical_independentVariants_4GenPC20ExpPC_2025-01-22/exported_matrix/coexpression_matrix_10497_genes.rds")

  gene_centers <- readRDS("gene_centers_16781_genes.rds")
  gene_deviations <- readRDS("gene_deviation_16781_genes.rds")

  message(sprintf("Maximum R2 in coexpression matrix: %.3f", max((coexp_10k^2)[upper.tri(coexp_10k, diag=F)])))
  adjacency_matrix <- (coexp_10k^2) > 0.01
  freq_high_coexp <- table(adjacency_matrix[upper.tri(adjacency_matrix, diag=F)])
  print(freq_high_coexp)

  # Only the diagonal is true. Code in if statement below is not executed.
  if (any(adjacency_matrix[upper.tri(adjacency_matrix, diag=F)])) {
    genes_after_pruning <- maximum_independent_set(adjacency_matrix, rownames(adjacency_matrix))
    print(length(genes_after_pruning))

    saveRDS(genes_after_pruning,sprintf("uncorrelated_genes_%s_r2max0.01.rds", length(genes_after_pruning)))

    coexp_10k_pruned <- coexp_10k[genes_after_pruning, genes_after_pruning]
    saveRDS(coexp_10k_pruned, sprintf("coexpression_matrix_%s_genes_pruned.rds", length(genes_after_pruning)))
    #fwrite(coexp_10k_pruned, sprintf("coexp_pca_%s_genes_prunes.rds", length(genes_after_pruning)))

    # Run EVD
    pca_on_gene_correlation_matrix_pruned <- pca_train(correlation_matrix = coexp_10k_pruned)

    saveRDS(pca_on_gene_correlation_matrix_pruned$eigenvectors, sprintf("coexp_pca_%s_genes_pruned_eigenvectors.rds", length(genes_after_pruning)))
    saveRDS(pca_on_gene_correlation_matrix_pruned$eigenvalues, sprintf("coexp_pca_%s_genes_pruned_eigenvalues.rds", length(genes_after_pruning)))

    write_pca_fit(eigenvectors, eigenvalues, gene_centers, gene_deviations, sprintf("coexp_pca_%s_genes_pruned", length(genes_after_pruning)))
  }

  # Run EVD
  pca_on_gene_correlation_matrix <- pca_train(correlation_matrix = coexp_10k)

  saveRDS(pca_on_gene_correlation_matrix$eigenvectors, "eigenvectors_10497_genes.rds")
  saveRDS(pca_on_gene_correlation_matrix$eigenvalues, "eigenvalues_10497_genes.rds")


  eigenvectors <- readRDS( "eigenvectors_10497_genes.rds")
  eigenvalues <- readRDS("eigenvalues_10497_genes.rds")

  gene_centers <- readRDS("gene_centers_16781_genes.rds")
  gene_deviations <- readRDS("gene_deviation_16781_genes.rds")

  out_label <- sprintf("coexp_pca_%s_genes", nrow(pca_on_gene_correlation_matrix$eigenvectors))
  dir.create(out_label)

  write_pca_fit(pca_on_gene_correlation_matrix$eigenvectors, pca_on_gene_correlation_matrix$eigenvalues,
                gene_centers, gene_deviations, file.path(out_label, out_label))


}

if (sys.nframe() == 0 && !interactive()) {
  main()
}