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
__program__ = "Annotate Loci"
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
class MafCalculator:
    def __init__(self, inclusion_path, maf_table, flipped,
                 table_name="filter_logs_full.log",
                 variant_inclusion_format="%s_SnpsToInclude.txt",
                 gene_inclusion_format="%s_GenesToInclude.txt"):
        self.overview_df = pd.read_table(os.path.join(inclusion_path, table_name), index_col=False)
        self.overview_df.set_index('Dataset', inplace=True)
        self.maf_table = maf_table[self.overview_df.index]
        self.maf_table[flipped] = 1 - self.maf_table.loc[flipped,:]
        self.overview_df['snp_inclusion_path'] = (
            self.overview_df.index.map(lambda name: os.path.join(inclusion_path, variant_inclusion_format % name)))
        self.overview_df['gene_inclusion_path'] = (
            self.overview_df.index.map(lambda name: os.path.join(inclusion_path, gene_inclusion_format % name)))
        self.snp_inclusion_df = self.load_inclusion_df('snp_inclusion_path')
        self.gene_inclusion_df = self.load_inclusion_df('gene_inclusion_path')
    def load_inclusion_df(self, column):
        # Create an empty dictionary to store dataframes
        dfs = {}
        # Generate dict of inclusion paths
        inclusion_paths = self.overview_df[column].to_dict()
        # Iterate over files in the directory
        for cohort, filepath in inclusion_paths.items():
            df = pd.read_csv(filepath)
            df[cohort] = 1
            # Set the index to the 'ID' column
            df.set_index('ID', inplace=True)
            # Add the dataframe to the dictionary
            dfs[cohort] = df
        # Merge the dataframes on their index
        merged_df = pd.concat(dfs.values(), axis=1)
        # Replace NaN values with 0 and convert to boolean
        merged_df = merged_df.fillna(0).astype(bool)
        return merged_df
    def calculate_maf(self, gene_variant_df):
        # Reformat the presence of the variants in the given dataframe
        variant_presence = (
            gene_variant_df
            .merge(self.snp_inclusion_df, right_index=True, left_on='variant', how='left')
            .set_index(['phenotype', 'variant']))
        print(variant_presence)
        # Reformat the presence of the genes in the given dataframe
        gene_presence = (
            gene_variant_df
            .merge(self.gene_inclusion_df, right_index=True, left_on='phenotype', how='left')
            .set_index(['phenotype', 'variant']))
        print(gene_presence)
        # Now that the tables displaying presence have both the same index, we can determine
        # for each combination if it is present or not.
        combined_presence = variant_presence & gene_presence
        print(combined_presence)
        # Now, reformat the maf table to also be according to this format.
        variant_maf = (
            gene_variant_df
            .merge(self.maf_table, right_index=True, left_on='variant', how='left')
            .set_index(['phenotype', 'variant']))
        # Now, for each cohort in the maf table, multiply all MAFs by the sample size
        variant_maf_weighted = variant_maf.mul(self.overview_df.loc[variant_maf.columns, "N"])
        print(variant_maf_weighted)
        # Now sum the weighted MAFs, and divide this by the total sample size
        maf = (variant_maf_weighted.where(combined_presence).sum(axis=1)
               / combined_presence.dot(self.overview_df["N"]))
        # Now return this
        return maf


class GencodeParser:
    def __init__(self, filename):
        self.filename = filename
        self.df = self.parse()
    def parse(self):
        genes = []
        with gzip.open(self.filename, 'rt') as f:
            for line in f:
                if not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if fields[2] == 'gene':
                        attributes = dict(item.split('=') for item in fields[8].split(';') if '=' in item)
                        genes.append({'gene_id': attributes['gene_id'],
                                      'chromosome': fields[0],
                                      'start': int(fields[3]),
                                      'end': int(fields[4])})
        df = pd.DataFrame(genes).astype({'start': 'Int64', 'end': 'Int64'})
        df['chromosome'] = df['chromosome'].str.lstrip('chr').replace({"M": "25", "X": "23", "Y": "24"}).astype('Int64')
        df['gene_id'] = df['gene_id'].str.split('.').str[0]
        print(df['chromosome'].unique())
        return df[['gene_id', 'chromosome', 'start', 'end']].drop_duplicates(keep='first')


# Functions

# Main
def main(argv=None):
    if argv is None:
        argv = sys.argv

    # Process input
    parser = argparse.ArgumentParser(description='Annotate significant eQTL results with variant and gene information')

    parser.add_argument('--input-file', dest='input_file', required=True,
                        help='Path to the table containing eQTL results')
    parser.add_argument('--variant-reference', dest='variant_reference', required=True,
                        help='Path to the table containing all SNPs from a reference panel')
    parser.add_argument('--gene-gff', dest='gene_gff', required=True,
                        help='Path to the Gencode GFF3 file')
    parser.add_argument('--maf-table', dest='maf', required=True,
                        help='Path to the table containing minor allele frequencies')
    parser.add_argument('--inclusion-path', dest='inclusion_path', required=True,
                        help='Inclusion_path')
    parser.add_argument('--out-prefix', dest='out_prefix', required=True,
                        help='Prefix to use for output file names')

    args = parser.parse_args(argv[1:])
    print(args)
    print("Loading variant reference from '{}'".format(args.variant_reference))

    variant_reference = (
        pd.read_csv(args.variant_reference, sep = ' ', dtype={'CHR': "Int64", 'bp': "Int64"})
        .drop(["allele1", "allele2"], axis=1)
        .rename({"ID": "variant", "bp": "bp_variant", "CHR": "chromosome_variant",
                 "str_allele1": "allele_ref", "str_allele2": "allele_eff"}, axis=1))

    print("Variant reference loaded:")
    print(variant_reference.head())

    print("Loading minor allele frequencies from '{}'".format(args.maf))

    maf_dataframe = (
        pd.read_table(args.maf)
        .drop(["MedianMaf", "CombinedMaf", "POS", "CHR"], axis=1)
        .rename({"ID": "variant", "OtherAllele": "other_allele_maf", "Allele": "allele_maf"}, axis=1)
        .set_index("variant").join(variant_reference))

    maf_dataframe = pd.merge(maf_dataframe, variant_reference,
             left_on="variant", right_on="variant", validate="1:1")

    maf_dataframe["flipped"] = maf_dataframe["allele_ref"] == maf_dataframe["allele_maf"]
    print((maf_dataframe["allele_ref"] == maf_dataframe["allele_maf"]).sum())
    print((maf_dataframe["allele_other"] == maf_dataframe["allele_maf"]).sum())
    assert np.alltrue(maf_dataframe["flipped"] == ~(maf_dataframe["allele_eff"] == maf_dataframe["allele_maf"]))

    print("Minor allele frequencies loaded:")
    print(maf_dataframe.head())

    print("Initiating MAF calculator...")

    maf_calculator = MafCalculator(
        inclusion_path=args.inclusion_path,
        maf_table=maf_dataframe,
        flipped=maf_dataframe["flipped"])

    print("MAF calculator initiated!")

    print("Loading gene annotations from '{}'".format(args.gene_gff))

    # gencode_parser = GencodeParser(args.gene_gff)
    # gene_dataframe = gencode_parser.df
    eqtls = pd.read_csv(args.input_file, sep = '\t')

    # Perform method
    # eqtls_annotated = (
    #     eqtls
    #     .merge(variant_reference, how="left", on="variant")
    #     .merge(gene_dataframe, how="left", left_on="phenotype", right_on="gene_id"))

    eqtls_annotated = (
        eqtls.merge(variant_reference, how="left", on="variant"))

    maf = maf_calculator.calculate_maf(eqtls_annotated[['variant', 'phenotype']])
    eqtls_annotated['allele_eff_freq'] = maf.values

    eqtls_annotated.to_csv("{}.csv.gz".format(args.out_prefix), sep="\t", index=False)
    # Output
    return 0


if __name__ == "__main__":
    sys.exit(main())
