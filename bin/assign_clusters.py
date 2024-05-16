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
def assign_clusters(input_file, max_size):
    current_cluster = []
    current_chrom = None
    current_start = None

    with open(input_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            chrom, start, end, name = fields[:4]
            start, end = int(start), int(end)

            if chrom != current_chrom or end - current_start > max_size:
                if current_cluster:
                    yield current_cluster
                current_cluster = [(chrom, start, end, name)]
                current_chrom = chrom
                current_start = start
            else:
                current_cluster.append((chrom, start, end, name))

    if current_cluster:
        yield current_cluster


# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    if len(argv) != 4:
        print("Usage: python script.py input.bed max_cluster_size output_prefix")
        sys.exit(1)

    input_file = sys.argv[1]
    max_size = int(sys.argv[2])
    output_prefix = sys.argv[3]

    for i, cluster in enumerate(assign_clusters(input_file, max_size), 1):
        with open("{}_cluster_{}.bed".format(output_prefix, i), 'w') as f:
            for region in cluster:
                f.write('\t'.join(map(str, [region[0], region[1], region[2], region[3], f"cluster_{i}"])) + '\n')


    return 0


if __name__ == "__main__":
    sys.exit(main())
