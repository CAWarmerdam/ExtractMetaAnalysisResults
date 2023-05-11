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


def find_uncorr_genes(r_squared_matrix, names, threshold=0.1):
    # Binarise the correlation matrix
    matrix_threshold = r_squared_matrix > threshold
    print(matrix_threshold)
    return maximum_independent_set(matrix_threshold, names)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Calculate correlations between z-scores and find uncorrelated genes')
    parser.add_argument('--zscores-file', dest='input_file', help='Path to the input CSV file',
                        default=None, required=False)
    parser.add_argument('--output-file', dest='output_file', help='Path to the output file')
    parser.add_argument('--gene-correlations', dest='gene_correlations', help='Path to gene correlations. Will ignore input_file.',
                        default=None, required=False)
    parser.add_argument('-t', '--threshold', type=float, default=0.5, help='Correlation threshold (default: 0.5)')
    args = parser.parse_args()

    corr_matrix = None

    if args.gene_correlations is not None:
        print("Loading gene correlations: {}".format(args.gene_correlations))
        corr_matrix = pd.read_csv(args.gene_correlations, index_col="phenotype")

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

        corr_matrix.to_csv("gene_correlation_matrix.csv.gz")

    if corr_matrix is None:
        raise ValueError("Correlation matrix has not been defined. "
                         "Provide either Z-score input file, or precalculated gene-gene correlations")

    abs_corr_matrix = corr_matrix.abs()

    print(abs_corr_matrix)

    uncorrelated_genes = find_uncorr_genes(abs_corr_matrix.to_numpy(),
                                           names=corr_matrix.columns, threshold=args.threshold)

    number_of_uncorrelated_genes = dict()

    for threshold in [x / 100.0 for x in range(5, 51, 5)]:
        print(threshold)
        uncorrelated_genes = find_uncorr_genes(abs_corr_matrix.to_numpy(),
                                               names=corr_matrix.columns, threshold=args.threshold)
        number_of_uncorrelated_genes[threshold] = len(uncorrelated_genes)

    pd.DataFrame.from_dict(number_of_uncorrelated_genes).to_csv("number_of_uncorrelated_genes.csv")

    # write the list of uncorrelated genes to a file
    with open(args.output_file, 'w') as f:
        f.write('\n'.join(sorted(uncorrelated_genes)))

    return 0


if __name__ == '__main__':
    sys.exit(main())
