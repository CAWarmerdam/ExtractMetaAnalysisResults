#!/usr/bin/env python3


import sys
import argparse
import pandas as pd


def maximum_independent_set(matrix):
    """
    Given a boolean matrix represented as a pandas dataframe, returns the vertices
    that are part of the maximum independent set using the algorithm that identifies
    the vertex with the least amount of edges, and then removes all vertices that have
    an edge with this vertex.
    """
    # Initialize the set of independent vertices to be empty
    independent_vertices = set()

    # Loop while the matrix is not empty
    while not matrix.empty:
        # Find the vertex with the least amount of edges
        least_connected = matrix.sum(axis=1).idxmin()
        print(least_connected)

        # Add the least connected vertex to the independent set
        independent_vertices.add(least_connected)

        # Remove all vertices that have an edge with the least connected vertex
        matrix = matrix.loc[matrix[least_connected] == 0, matrix[least_connected] == 0]
        print(matrix.shape)

    return independent_vertices


def find_uncorr_genes(r_squared_matrix, threshold):
    # Binarise the correlation matrix
    matrix_threshold = pd.DataFrame(
        (r_squared_matrix >= threshold),
        index=r_squared_matrix.index,
        columns=r_squared_matrix.columns)
    print(matrix_threshold)
    return maximum_independent_set(matrix_threshold)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Calculate correlations between z-scores and find uncorrelated genes')
    parser.add_argument('input_file', help='Path to the input CSV file')
    parser.add_argument('output_file', help='Path to the output file')
    parser.add_argument('-t', '--threshold', type=float, default=0.5, help='Correlation threshold (default: 0.5)')
    args = parser.parse_args()

    # read the input file into a pandas DataFrame
    df = pd.read_csv(args.input_file, sep="\t")

    print(df)

    # pivot the DataFrame to get a matrix of z-scores
    matrix = df.pivot(index='variant', columns='phenotype', values='z_score')

    print(matrix)

    # calculate the pairwise correlations between genes
    corr_matrix = matrix.corr()

    corr_matrix.to_csv("gene_correlation_matrix.csv.gz")

    r_squared_matrix = corr_matrix.pow(2)

    print(r_squared_matrix)

    uncorrelated_genes = find_uncorr_genes(r_squared_matrix, args.threshold)

    # write the list of uncorrelated genes to a file
    with open(args.output_file, 'w') as f:
        f.write('\n'.join(sorted(uncorrelated_genes)))

    return 0


if __name__ == '__main__':
    sys.exit(main())
