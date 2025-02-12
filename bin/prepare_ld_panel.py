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
import os
import re
import sys
import glob

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import stats

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


# Constants
PYARROW_SCHEMA_LD = pa.schema([
    ("phenotype", pa.string()),
    ("variant_index", pa.int64()),
    ("rho", pa.float64())])


# Classes

# Functions
def load_genes(file_path):
    """Load the list of genes from a file."""
    with open(file_path, 'r') as f:
        return [line.strip() for line in f]


def t_to_cor(t_stats, df):
    # Compute the correlation coefficients
    r = t_stats / np.sqrt(df + t_stats**2)
    return r


def get_variant_indices(reference_file, chromosome, min, max):
    """Get the min and max variant indices for a given chromosome."""
    variant_reference = pq.read_table(reference_file)
    chromosome_data = variant_reference.filter(pa.compute.equal(variant_reference.column('chromosome'), chromosome)).to_pandas()
    return chromosome_data[np.logical_and(chromosome_data.variant_index >= min, chromosome_data.variant_index <= max)]


def process_and_write_gene_data(dataset_folder, gene, min_index, max_index):
    """Filter data for a gene, perform calculations, and write results to Parquet."""
    gene_file_paths = glob.glob(os.path.join(dataset_folder, f"phenotype={gene}/*.parquet"))
    assert len(gene_file_paths) == 1
    gene_file_path = gene_file_paths[0]
    if not os.path.exists(gene_file_path):
        print(f"Warning: File for gene {gene} not found. Skipping.")
        return
    filters = [[('variant_index', '>=', min_index), ('variant_index', '<=', max_index)]]
    # Load parquet data for gene
    summary_stats = pq.ParquetDataset(gene_file_path, validate_schema=True, filters=filters).read().to_pandas()
    # Calculate rho
    t_stat = summary_stats['beta'] / summary_stats['standard_error']
    df = summary_stats['sample_size'] - 1
    rho = t_to_cor(t_stat, df)
    return summary_stats['variant_index'], rho


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv
    # Process input
    parser = argparse.ArgumentParser(description="Summary statistics to prepare the LD dataset.")
    parser.add_argument("--dataset-folder", required=True, help="Path to the partitioned parquet dataset folder.")
    parser.add_argument("--genes-file", required=True, help="Path to the file containing the list of genes.")
    parser.add_argument("--variant-reference", required=True, help="Path to the variant reference parquet file.")
    parser.add_argument("--output-folder", required=True, help="Output folder.")
    parser.add_argument("--chromosome", type=int, required=True, help="Chromosome number to process.")
    parser.add_argument("--min-max", type=int, nargs=2, required=True, help="Variant index range to process")
    args = parser.parse_args()

    genes = load_genes(args.genes_file)[:20]
    n_genes = len(genes)

    k_variants = 500000

    max_variant_index = args.min_max[1]
    min_variant_index = args.min_max[0]
    chr_output_folder = os.path.join(args.output_folder, f"chr={args.chromosome}")
    output_file = os.path.join(chr_output_folder, f"{min_variant_index}_{max_variant_index}.parquet")

    chr_variant_indices = get_variant_indices(args.variant_reference, args.chromosome, min_variant_index, max_variant_index)
    chr_chunks = np.array_split(chr_variant_indices, np.ceil(chr_variant_indices.shape[0] / k_variants))

    os.makedirs(chr_output_folder, exist_ok=True)

    pyarrow_ld_schema = pa.schema([("variant_index", pa.int64())] + [(gene_id, pa.float64()) for gene_id in genes])
    parquet_writer = pq.ParquetWriter(output_file, pyarrow_ld_schema)
    print(f"Nr of chunks: {len(chr_chunks)}")

    for chr_chunk in chr_chunks:
        result_table = None

        min_index = chr_chunk.variant_index.min()
        max_index = chr_chunk.variant_index.max()

        i = 0
        for gene in genes:
            variant_indices, rho = process_and_write_gene_data(args.dataset_folder, gene, min_index, max_index)

            if not result_table:
                result_table = pa.Table.from_pydict({
                    'variant_index': variant_indices,
                    gene: rho
                })
            else:
                assert np.array_equal(result_table.column('variant_index').to_numpy(), variant_indices)
                result_table = result_table.append_column(gene, pa.array(rho))
            n_variants = len(variant_indices)
            i+=1
            print(f"Chunk: {min_index}-{max_index}, {i}/{n_genes} done, found {n_variants} variants", end="\r")

        print(f"Chunk: {min_index}-{max_index}, {i}/{n_genes} done, found {n_variants} variants")
        print("Writing results")
        parquet_writer.write_table(result_table)
        print("Chunk Done!")

    parquet_writer.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
