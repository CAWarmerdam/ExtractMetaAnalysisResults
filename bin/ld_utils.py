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


class LdCalculator:
    def __init__(self, parquet_dataset, variant_reference, observations=None):
        self.reference = variant_reference
        self.parquet_dataset = parquet_dataset
        self.observations = observations

    def _calc_ld(self, x):
        """
        Extracts a locus as a correlation matrix (LD matrix) from a Parquet dataset.

        Args:
            x (pd.DataFrame): The input dataset containing 'variant_index', and values across columns
            variant_index_start (int): The starting index for filtering variants.
            variant_index_end (int): The ending index for filtering variants.

        Returns:
            np.ndarray: The computed LD matrix.
        """
        rho_mat = x.to_numpy()

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

        return ld_matrix

    def calculate_ld_non_continuous(self, locus_bed):
        start_time = time.time()

        locus_grouped = (locus_bed.groupby((locus_bed.end.shift() - locus_bed.start).lt(0).cumsum())
                         .agg({'chromosome': 'first', 'start': 'first', 'end': 'last', 'name': 'sum'}))
        print("Starting to calculate LD:")
        print(locus_grouped)

        x_partial_list = list()

        for index, continuous_locus in locus_grouped.iterrows():
            # Extract locus information
            locus_chromosome = continuous_locus["chromosome"].unique()[0]
            locus_start = continuous_locus["start"].min()
            locus_end = continuous_locus["end"].max()
            print(f"    - Extracting input for LD calculations: {locus_chromosome}:{locus_start}-{locus_end}")

            # Get the variants in the locus
            locus_variant_reference = self.reference[
                (self.reference["chromosome"] == locus_chromosome) &
                (self.reference["bp"].between(locus_start, locus_end))
                ]

            # Get the indices of the variants
            variant_index_start = locus_variant_reference["variant_index"].min()
            variant_index_end = locus_variant_reference["variant_index"].max()

            # Filter dataset based on variant index range
            x_loc = self.get_x_for_locus(variant_index_end, variant_index_start)
            x_partial_list.append(x_loc)
            print(f"    - Found {x_loc.shape[0]} variants")

        x = pd.concat(x_partial_list, axis=0)
        print(f"Found {x.shape[0]} total variants!")
        ld_matrix = self._calc_ld(x)
        end_time = time.time()
        print("LD calculation done!")
        print(f"Time taken: {end_time - start_time} seconds")
        return ld_matrix, x.index.values

    def get_x_for_locus(self, variant_index_end, variant_index_start):
        if self.observations is not None:
            filters = [[
                ("phenotype", 'in', self.observations),
                ("variant_index", ">=", variant_index_start),
                ("variant_index", "<=", variant_index_end)]]
        else:
            filters = [[
                ("variant_index", ">=", variant_index_start),
                ("variant_index", "<=", variant_index_end)]]
        return (pq.ParquetDataset(self.parquet_dataset,
                                  filters=filters).read().to_pandas()
                .pivot(index="variant_index", columns="phenotype", values="rho").fillna(0))


class LdCalculatorWide(LdCalculator):
    def get_x_for_locus(self, variant_index_end, variant_index_start):
        if self.observations is not None:
            filters = [[
                ("phenotype", 'in', self.observations),
                ("variant_index", ">=", variant_index_start),
                ("variant_index", "<=", variant_index_end)]]
        else:
            filters = [[
                ("variant_index", ">=", variant_index_start),
                ("variant_index", "<=", variant_index_end)]]
        return pq.ParquetDataset(self.parquet_dataset, filters=filters).read().to_pandas()



def main(argv=None):
    if argv is None:
        argv = sys.argv
    return 0


if __name__ == "__main__":
    sys.exit(main())
