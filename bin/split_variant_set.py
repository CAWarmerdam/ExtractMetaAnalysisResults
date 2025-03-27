#!/usr/bin/env python3

"""
Created:      04/12/2024
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

# Standard imports.
import argparse
import sys

import numpy as np
import pandas as pd

import pyarrow.parquet as pq


# Metadata
__program__ = "CNV-caller"
__author__ = "C.A. (Robert) Warmerdam"
__email__ = "c.a.warmerdam@umcg.nl"
__license__ = "GPLv3"
__version__ = 1.0
__description__ = "{} is a program developed and maintained by {}. " \
                  "This program is licensed under the {} license and is " \
                  "provided 'as-is' without any warranty or indemnification " \
                  "of any kind.".format(__program__,
                                        __author__,
                                        __license__)


# Classes

# Functions
def split_dataframe_by_group(df, group_col, k):
    variants_per_chunk = np.ceil(df.shape[0] / k)

    # Print variants_per_chunk
    print(f"Using k={k}, n={df.shape[0]}")
    print(f"Variants per chunk: {variants_per_chunk}")

    # Group by the specified column
    grouped = list(df.groupby(group_col))

    k_left = k

    chunks = {
        "chromosome": list(),
        "min_variant_index": list(),
        "max_variant_index": list(),
        "chunk_id": list()
    }

    # Distribute groups to chunks in a round-robin fashion (or another strategy)
    for chr_name, chr_df in grouped:
        # Find the chunk with the least rows
        k_chromosome = np.round(chr_df.shape[0] / variants_per_chunk)
        print(f"Splitting chromosome {chr_name}: {chr_df.shape[0]} rows")
        print(f"Taking minimum of remaining chunks ({k_left}), and {k_chromosome}")
        chr_chunks = np.array_split(
            chr_df,
            min(k_chromosome, k_left))
        print(f"Splitted in to {len(chr_chunks)} chunks")
        k_left -= len(chr_chunks)
        print(f"{k_left}/{k} remaining.")
        for i, chunk in enumerate(chr_chunks):
            if chunk.empty:
                continue

            chunks["chromosome"].append(chr_name)
            chunks["chunk_id"].append(f"{chr_name}_{i}")
            chunks["min_variant_index"].append(chunk["variant_index"].min())
            chunks["max_variant_index"].append(chunk["variant_index"].max())

    # Convert list of dicts to output dataframe
    return pd.DataFrame.from_dict(chunks)

# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv
    # Process input
    parser = argparse.ArgumentParser(description="Split variant reference into chunks.")
    parser.add_argument("--variant-reference", required=True, help="Path to the variant reference parquet file.")
    parser.add_argument("--chunks", required=True, type=int, help="Number of chunks")
    parser.add_argument("--output", required=True, help="Output file")
    args = parser.parse_args()

    variant_reference = pq.read_table(args.variant_reference).to_pandas()

    chunks = split_dataframe_by_group(variant_reference, 'chromosome', args.chunks)

    chunks.to_csv(args.output, index=False)
    return 0


if __name__ == "__main__":
    sys.exit(main())
