#!/usr/bin/env python3

"""
Created:      29/11/2024
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2024 C.A. Warmerdam

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""


import sys
import pandas as pd
import pyarrow.parquet as pq
import numpy as np
import time


def get_ld_matrix(permuted_dataset_path, phenotype_list, variant_index_start, variant_index_end):
    """
    Extracts a locus as a correlation matrix (LD matrix) from a Parquet dataset.

    Args:
        permuted_dataset (pyarrow): The input dataset containing 'variant_index', 'phenotype', and 'rho' columns.
        phenotype_list (list):List of phenotypes to filter on.
        variant_index_start (int): The starting index for filtering variants.
        variant_index_end (int): The ending index for filtering variants.

    Returns:
        np.ndarray: The computed LD matrix.
    """
    # Filter dataset based on variant index range
    df_filtered = (pq.ParquetDataset(
        permuted_dataset_path,
        filters=[[
            ("phenotype", 'in', phenotype_list),
            ("variant_index", ">=", variant_index_start),
            ("variant_index", "<=", variant_index_end)]]).read().to_pandas()
                   .pivot(index="variant_index", columns="phenotype", values="rho").fillna(0))

    variant_names = df_filtered.index.values

    rho_mat = df_filtered.to_numpy()

    start_time = time.time()

    # Centering the matrix (subtract row means)
    rho_mat -= rho_mat.mean(axis=1, keepdims=True)

    # Standardizing each row
    row_norms = np.sqrt(np.sum(rho_mat**2, axis=1, keepdims=True))
    rho_mat /= row_norms

    # Compute LD matrix
    ld_matrix = np.dot(rho_mat, rho_mat.T)

    end_time = time.time()
    print(f"Time taken: {end_time - start_time:.4f} seconds")

    return ld_matrix, variant_names


def main(argv=None):
    if argv is None:
        argv = sys.argv
    return 0


if __name__ == "__main__":
    sys.exit(main())
