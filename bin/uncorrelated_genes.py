#!/usr/bin/env python3


import sys
import argparse
import pandas as pd
import numpy as np
import numpy.ma as ma


def maximum_independent_set(dataframe):
    """
    Given a boolean matrix represented as a pandas dataframe, returns the vertices
    that are part of the maximum independent set using the algorithm that identifies
    the vertex with the least amount of edges, and then removes all vertices that have
    an edge with this vertex.
    """
    # Initialize the set of independent vertices to be empty
    independent_vertices = set()

    # Convert pandas DataFrame to numpy array
    matrix = dataframe.to_numpy()

    # Sample names
    names = dataframe.columns.to_numpy()

    # Loop while the matrix is not empty
    while not matrix.size == 0:
        # Find the vertex with the least amount of edges
        least_connected = matrix.sum(axis=1).argmin()
        print(names[least_connected])

        # Add the least connected vertex to the independent set
        independent_vertices.add(names[least_connected])

        # Remove all vertices that have an edge with the least connected vertex
        no_edges = matrix[least_connected] == 0
        matrix = matrix[no_edges, no_edges]
        names = names[no_edges]
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
    corr_matrix = ma.corrcoef(ma.masked_invalid(matrix.to_numpy()), rowvar=False)

    # Make a pandas dataframe
    corr_dataframe = pd.DataFrame(corr_matrix, index=matrix.columns, columns=matrix.columns)

    corr_dataframe.to_csv("gene_correlation_matrix.csv.gz")

    r_squared_matrix = corr_dataframe.pow(2)

    print(r_squared_matrix)

    uncorrelated_genes = find_uncorr_genes(r_squared_matrix, args.threshold)

    # write the list of uncorrelated genes to a file
    with open(args.output_file, 'w') as f:
        f.write('\n'.join(sorted(uncorrelated_genes)))

    return 0


if __name__ == '__main__':
    sys.exit(main())
