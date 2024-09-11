#!/usr/bin/env Rscript


# Load libraries
library(argparse)
library(tidyverse)
library(data.table)
library(RSpectra)
library(rhdf5)

# Declare constants

# Declare function definitions

ZtoCorLude<- function(z, n) {
  mathLog10N <- log10(n)
  corr <- z * 10^(0.00332727459306132 - 0.500803673176124 * mathLog10N) -
    z^3 * 10^(-0.600314180613649 - 1.50039529106925 * mathLog10N) +
    z^5 * 10^(-1.3415016045602 - 2.48579979106749 * mathLog10N) -
    z^7 * 10^(-2.53051421275022 - 3.39996348904233 * mathLog10N)
  return(corr)
}

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
  logp <- min(pt(t, df, log.p = T), pt(t, df, lower.tail=FALSE, log.p = T))#normaly need *2 since this is two-side. However only used to covert to z-score
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
  parser$add_argument('--input-prefix', dest='input_prefix', required=F)
  parser$add_argument('--output-file', dest='output_file', help='Path to the output file')
  parser$add_argument('--gene-correlations', dest='gene_correlations', help='Path to gene correlations. Will ignore input_file.',
                      default=NULL, required=F)
  parser$add_argument('-c', '--n-components', type="double", default=1000, help='Number of components to write')
  parser$add_argument('-n', '--n-threshold', type="double", default=NULL, help='Minimum sample size of gene-variant pairs to consider for gene-gene correlations')
  args <- parser$parse_args()

  corr_matrix <- NULL
  sample_size_threshold <- args$n_threshold

  if (!is.null(args$gene_correlations)) {
    cat("Loading gene correlations:", args$gene_correlations, "\n")
    corr_matrix <- read.csv(args$gene_correlations, row.names = 1, sep = "\t")
  }

  input_prefix <- args$input_prefix
  if (!is.null(input_prefix)) {
    if (!is.null(corr_matrix)) {
      stop("Correlation matrix is already defined, skipping input prefix")
    }

    input_file_z <- paste0(input_prefix, ".z.txt.gz")
    input_file_n <- paste0(input_prefix, ".n.txt.gz")
    cat("Loading Z-scores:", input_file_z, "\n")
    z_matrix <- as.matrix(fread(input_file_z, sep = "\t"), rownames = 1)

    cat("Loading sample-sizes:", input_file_n, "\n")
    n_matrix <- as.matrix(fread(input_file_n, sep = "\t"), rownames = 1)

    if (!is.null(sample_size_threshold)) {
      variants_where_all_na <- rowSums(n_matrix < max(n_matrix)) == ncol(z_matrix)
    } else {
      variants_where_all_na <- rowSums(n_matrix < sample_size_threshold) == ncol(z_matrix)
    }

    all(colnames(n_matrix) == colnames(z_matrix))
    all(rownames(n_matrix) == rownames(z_matrix))

    n_matrix <- n_matrix[!variants_where_all_na,]
    z_matrix <- z_matrix[!variants_where_all_na,]

    # genes_where_all_na <- colSums(is.na(z_matrix)) == nrow(z_matrix)
    #
    # n_matrix <- n_matrix[,!genes_where_all_na]
    # z_matrix <- z_matrix[,!genes_where_all_na]
    #
    # variants_without_na <- rowSums(is.na(z_matrix)) == 0
    # n_matrix <- n_matrix[variants_without_na,]
    # z_matrix <- z_matrix[variants_without_na,]

    z_as_corr_matrix <- t(zToCor(z_matrix, n_matrix-1))

    # Center each variable
    z_as_corr_matrix <- z_as_corr_matrix - rowMeans(z_as_corr_matrix)
    # Standardize each variable
    z_as_corr_matrix <- z_as_corr_matrix / sqrt(rowSums(z_as_corr_matrix^2))
    # Calculate correlations
    corr_matrix <- tcrossprod(z_as_corr_matrix)

    fwrite(corr_matrix, "gene_correlation_matrix.csv.gz", sep = "\t", col.names = T, row.names = T, quote = FALSE, na = "NA")
    saveRDS(corr_matrix, file = "gene_correlation_matrix.rds")
    message("Done!")
  }

  if (is.null(corr_matrix)) {
    stop("Correlation matrix has not been defined. Provide either Z-score input file, or precalculated gene-gene correlations")
  }

  # Perform eigen decomposition
  eigendecomposition <- RSpectra::eigs_sym(corr_matrix, k=1000)
  # Extract the eigenvectors from the eigen decomposition
  eigenvectors <- eigendecomposition$vectors
  # Assign rownames to the eigenvectors
  rownames(eigenvectors) <- rownames(corr_matrix)
  # Save the eigenvectors as an RDS file
  fwrite(eigenvectors, "gene_eigenvectors.csv.gz", sep = "\t", col.names = T, row.names = T, quote = FALSE, na = "NA")
  saveRDS(eigenvectors, file = "gene_eigenvectors.rds")

  # Process output
}

if (sys.nframe() == 0 && !interactive()) {
  main()
}