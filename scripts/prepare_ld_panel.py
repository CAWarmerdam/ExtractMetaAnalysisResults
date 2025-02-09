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


def get_variant_indices(reference_file, chromosome):
    """Get the min and max variant indices for a given chromosome."""
    variant_reference = pq.read_table(reference_file)
    chromosome_data = variant_reference.filter(pa.compute.equal(variant_reference.column('chromosome'), chromosome))
    return chromosome_data


def process_and_write_gene_data(dataset_folder, gene, min_index, max_index, parquet_writer):
    """Filter data for a gene, perform calculations, and write results to Parquet."""
    gene_file_path = os.path.join(dataset_folder, f"phenotype_{gene}.parquet")
    if not os.path.exists(gene_file_path):
        print(f"Warning: File for gene {gene} not found. Skipping.")
        return
    filters = [('variant_index', '>=', min_index), ('variant_index', '<=', max_index)]
    # Load parquet data for gene
    summary_stats = pq.ParquetDataset(gene_file_path, filters=filters).read().to_pandas()
    # Calculate rho
    t_stat = summary_stats['beta'] / summary_stats['standard_error']
    df = summary_stats['sample_size'] - 1
    rho = t_to_cor(t_stat, df)
    # Generate new output table
    result_table = pa.Table.from_pydict({
        'gene': [gene] * summary_stats.shape[0],
        'variant_index': summary_stats['variant_index'],
        'calculation_result': rho
    })

    # Write results to the output folder
    parquet_writer.write_table(result_table)


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
    args = parser.parse_args()

    genes = load_genes(args.genes_file)

    k_variants = 10000

    chr_variant_indices = get_variant_indices(args.variant_reference, args.chromosome)
    chr_chunks = np.array_split(chr_variant_indices, np.ceil(chr_variant_indices.shape[0] / k_variants))

    output_file = os.path.join(args.output_folder, f"ld_panel_chr{args.chromosome}.parquet")

    os.makedirs(args.output_folder, exist_ok=True)
    parquet_writer = pq.ParquetWriter(output_file, PYARROW_SCHEMA_LD)

    for gene in genes:
        process_and_write_gene_data(args.dataset_folder, gene, min_index, max_index, parquet_writer)

    parquet_writer.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
