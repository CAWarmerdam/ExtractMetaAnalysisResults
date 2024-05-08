#!/usr/bin/env python3


import sys
import argparse
import pandas as pd
import numpy as np
import numpy.ma as ma


def maximum_independent_set(matrix, names):
    """
    Given a boolean matrix represented as a pandas dataframe, returns the vertices
    that are part of the maximum independent set using the algorithm that identifies
    the vertex with the least amount of edges, and then removes all vertices that have
    an edge with this vertex.
    """
    # Initialize the set of independent vertices to be empty
    independent_vertices = set()

    # First collect every gene that is totally independent.
    # Removing these will not have an effect on other genes.
    unconnected_vertices = matrix.sum(axis=1) == 1
    independent_vertices.update(set(names[unconnected_vertices]))

    # Now remove these from the matrix and the names
    matrix = matrix[np.ix_(~unconnected_vertices, ~unconnected_vertices)]
    names = names[~unconnected_vertices]

    # Loop while the matrix is not empty
    while not matrix.size == 0:
        # Find the vertex with the least amount of edges
        least_connected = matrix.sum(axis=1).argmin()
        print(names[least_connected])

        # Add the least connected vertex to the independent set
        independent_vertices.add(names[least_connected])

        # Remove all vertices that have an edge with the least connected vertex
        no_edges = matrix[least_connected, :] == 0

        matrix = matrix[np.ix_(no_edges, no_edges)]
        names = names[no_edges]
        print(matrix.shape)

    return independent_vertices


def find_uncorr_genes(r_abs_matrix, names, threshold=0.1):
    # Binarise the correlation matrix
    matrix_threshold = r_abs_matrix > threshold
    print(matrix_threshold)
    return maximum_independent_set(matrix_threshold, names)


def write_uncorrelated_genes(uncorrelated_genes, output_file):
    # write the list of uncorrelated genes to a file
    with open(output_file, 'w') as f:
        f.write('\n'.join(sorted(uncorrelated_genes)))


def grid_search(abs_corr_matrix):
    number_of_uncorrelated_genes = dict()
    print("Threshold\tN", file=open('number_of_uncorrelated_genes.csv', 'w'))

    for threshold in [x / 100.0 for x in range(5, 20, 1)]:
        print(threshold)
        uncorrelated_genes = find_uncorr_genes(abs_corr_matrix.to_numpy(),
                                               names=corr_matrix.columns, threshold=threshold)
        number_of_uncorrelated_genes[threshold] = len(uncorrelated_genes)
        print("{}\t{}".format(threshold, len(uncorrelated_genes)), file=open('number_of_uncorrelated_genes.csv', 'a'))

        # write the list of uncorrelated genes to a file
        with open("uncorrelated_genes_N{}.txt".format(len(uncorrelated_genes)), 'w') as f:
            f.write('\n'.join(uncorrelated_genes))


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Calculate correlations between z-scores and find uncorrelated genes')
    parser.add_argument('--zscores-file-long', dest='input_file', help='Path to the input CSV file',
                        default=None, required=False)
    parser.add_argument('--input-prefix', dest='input_prefix', required=False)
    parser.add_argument('--output-file', dest='output_file', help='Path to the output file')
    parser.add_argument('--gene-correlations', dest='gene_correlations', help='Path to gene correlations. Will ignore input_file.',
                        default=None, required=False)
    parser.add_argument('-t', '--threshold', type=float, default=0.5, help='Correlation threshold (default: 0.5)')
    parser.add_argument('-n', '--n-threshold', type=float, default=0, help='Minimum sample size of gene-variant pairs to consider for gene-gene correlations')
    args = parser.parse_args()

    corr_matrix = None
    sample_size_threshold = args.n_threshold

    if args.gene_correlations is not None:
        print("Loading gene correlations: {}".format(args.gene_correlations))
        corr_matrix = pd.read_csv(args.gene_correlations, index_col=0, sep="\t")

    if args.input_file is not None:
        if corr_matrix is not None:
            raise ValueError("Correlation matrix is already defined, skipping Z-score input file")

        print("Loading Z-scores: {}".format(args.input_file))
        # read the input file into a pandas DataFrame
        df = pd.read_csv(args.input_file, sep="\t")

        print(df)

        df.drop_duplicates(inplace=True)

        # pivot the DataFrame to get a matrix of z-scores
        matrix = df.pivot(index='variant', columns='phenotype', values='z_score')

        print(matrix)

        # calculate the pairwise correlations between genes
        corr_matrix = matrix.corr()

        corr_matrix.to_csv("gene_correlation_matrix.csv.gz", sep="\t", index_label="phenotype")

    if args.input_prefix is not None:
        if corr_matrix is not None:
            raise ValueError("Correlation matrix is already defined, skipping input prefix")

        input_file_z = args.input_prefix + ".z.txt"
        input_file_n = args.input_prefix + ".n.txt"
        print("Loading Z-scores: {}".format(input_file_z))
        # read the input file into a pandas DataFrame
        z_matrix = pd.read_csv(input_file_z, sep="\t", index_col='variant')
        print("Loading sample-size threshold: {}".format(input_file_n))
        n_matrix = pd.read_csv(input_file_n, sep="\t", index_col='variant')

        print(z_matrix)

        # calculate the pairwise correlations between genes
        filtered_z_matrix = (
            z_matrix[n_matrix >= sample_size_threshold]
            .dropna(axis=0, how='all', inplace=False))
        corr_matrix = filtered_z_matrix.loc[:,filtered_z_matrix.isna().sum() < len(filtered_z_matrix) * 0.5].corr()

        corr_matrix.to_csv("gene_correlation_matrix.csv.gz", sep="\t", index_label="phenotype")

    if corr_matrix is None:
        raise ValueError("Correlation matrix has not been defined. "
                         "Provide either Z-score input file, or precalculated gene-gene correlations")

    output_file = args.output_file
    threshold = args.threshold

    abs_corr_matrix = corr_matrix.abs()

    print(abs_corr_matrix)

    uncorrelated_genes = find_uncorr_genes(abs_corr_matrix.to_numpy(),
                                           names=corr_matrix.columns, threshold=threshold)

    write_uncorrelated_genes(uncorrelated_genes, output_file)

    grid_search(abs_corr_matrix)
    print("Done")

    return 0


if __name__ == '__main__':
    sys.exit(main())
