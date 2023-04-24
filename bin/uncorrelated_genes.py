#!/usr/bin/env python2

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

        # Add the least connected vertex to the independent set
        independent_vertices.add(least_connected)

        # Remove all vertices that have an edge with the least connected vertex
        matrix = matrix.loc[matrix[least_connected] == 0, matrix[least_connected] == 0]

    return independent_vertices


def find_uncorr_genes(corr_matrix, threshold):
    # Binarise the correlation matrix
    return maximum_independent_set(corr_matrix.abs() >= threshold)


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
    df = pd.read_csv(args.input_file)

    # pivot the DataFrame to get a matrix of z-scores
    matrix = df.pivot(index='phenotype', columns='variant', values='z_score')

    # calculate the pairwise correlations between genes
    corr_matrix = matrix.corr()

    uncorrelated_genes = find_uncorr_genes(corr_matrix, args.threshold)

    # write the list of uncorrelated genes to a file
    with open(args.output_file, 'w') as f:
        f.write('\n'.join(sorted(uncorrelated_genes)))

    return 0


if __name__ == '__main__':
    main(sys.exit())