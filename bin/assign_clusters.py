#!/usr/bin/env python3

"""
Created:      14/02/2024
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
import os
import sys
import argparse

import pandas as pd

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

# Classes

# Functions
def assign_clusters(bed, max_size, max_n):
    current_cluster = []
    current_chrom = None
    current_start = None

    for name, line in bed.groupby(['name', 'gene_clusters'], sort=False):
        chrom = line["chromosome"].values[0]
        start = line["start"].min()
        end = line["end"].max()
        name = line["name"].values[0]
        start, end = int(start), int(end)
        print(f"{chrom}:{start}-{end}: {name}, {len(line['name'])} nearby loci")
        print(line)

        if current_start is not None and (chrom != current_chrom or (end - current_start) > max_size or len(current_cluster) >= max_n):
            print(chrom != current_chrom)
            print((end - current_start) > max_size)
            print(len(current_cluster) >= max_n)
            print(f"Finalizing cluster")
            if current_cluster:
                yield pd.concat(current_cluster, axis=0)
            print(f"Wrote cluster of {len(current_cluster)} genes")
            print(current_cluster)
            print("---")
            current_cluster = [line]
            current_chrom = chrom
            current_start = start
        else:
            current_cluster.append(line)
            if current_chrom is None and current_start is None:
                current_chrom = chrom
                current_start = start

    if current_cluster:
        yield pd.concat(current_cluster, axis=0)


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description="Assign clusters to bed file based on maximum distances and a max number of genes")

    parser.add_argument("--input-file", type=str, help="Path to the input bed file.")
    parser.add_argument("--max-distance", type=int, help="Maximum distance between first gene in cluster and last gene of cluster (including their gene bodies; from the TSS of gene 1 to the TES of gene n).")
    parser.add_argument("--output-prefix", type=str, help="Prefix for the output files")
    parser.add_argument("--max-n", type=int, help="Maximum number of items")

    args = parser.parse_args(argv[1:])

    input_file = args.input_file
    max_d = args.max_distance
    output_prefix = args.output_prefix
    max_n = args.max_n

    print(f"Input File: {input_file}")
    print(f"Max Size: {max_d}")
    print(f"Output Prefix: {output_prefix}")
    print(f"Max N: {max_n}")

    bed = pd.read_table(input_file, header=None, names=["chromosome", "start", "end", "name", "gene_clusters"])
    bed.sort_values('start', inplace=True)

    cluster_df = assign_clusters(bed, max_d, max_n)

    for i, cluster in enumerate(cluster_df, 1):
        cluster["locus"] = i
        cluster.to_csv("{}_cluster_{}.bed".format(output_prefix, i), sep='\t', header=False, index=False)

    return 0


if __name__ == "__main__":
    sys.exit(main())
