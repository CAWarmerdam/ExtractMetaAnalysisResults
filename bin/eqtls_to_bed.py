#!/usr/bin/env python3

"""
Created:      15/04/2023
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2023 C.A. Warmerdam

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
import gzip

import numpy as np
import pandas as pd

# Metadata
__program__ = "eQTLs2Bed"
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

# Main
from annotate_loci import GencodeParser


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    parser = argparse.ArgumentParser(description='Annotate significant eQTL results with variant and gene information')

    parser.add_argument('--input-file', dest='input_file', required=True,
                        help='Path to the table containing significant eQTL results')
    parser.add_argument('--variant-reference', dest='variant_reference', required=True,
                        help='Path to the table containing all SNPs from a reference panel')
    parser.add_argument('--gene-ggf', dest='gene_ggf', required=True,
                        help='Path to the Gencode GFF3 file')
    parser.add_argument('--out-prefix', dest='out_prefix', required=True,
                        help='Prefix to use for output file names')

    args = parser.parse_args(argv[1:])
    print(args)

    variant_reference = (
        pd.read_csv(args.variant_reference, sep = ' ', dtype={'CHR': "Int64", 'bp': "Int64"})
        .drop(["allele1", "allele2", "str_allele1", "str_allele2"], axis=1)
        .rename({"ID": "variant", "bp": "bp_variant", "CHR": "chromosome_variant"}, axis=1))
    gencode_parser = GencodeParser(args.gene_ggf)
    gene_dataframe = gencode_parser.df
    eqtls = pd.read_csv(args.input_file, sep = '\t')

    # Perform method
    eqtls_annotated = (
        eqtls
        .merge(variant_reference, how="left", on="variant")
        .merge(gene_dataframe, how="left", left_on="phenotype", right_on="gene_id"))

    # Identify genes that have a cis-effect
    preselection_cis = \
        eqtls_annotated.loc[(eqtls_annotated.chromosome_variant == eqtls_annotated.chromosome), :]

    # Select genes for which we can find a cis-effect
    cis_genes = preselection_cis.loc[np.logical_or(
        (preselection_cis.bp_variant - preselection_cis.start).abs() < 1*10**6,
        (preselection_cis.bp_variant - preselection_cis.end).abs() < 1*10**6),
                                     ["chromosome", "start", "end", "phenotype"]].drop_duplicates()

    # Select all significant variants
    variants = eqtls_annotated.loc[:,["chromosome_variant", "bp_variant", "bp_variant", "phenotype"]]

    # Output bed files for which to
    cis_genes.to_csv(".".join([args.out_prefix, "genes.bed"]), sep='\t', header=False, index=False)
    variants.to_csv(".".join([args.out_prefix, "variants.bed"]), sep='\t', header=False, index=False)

    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
