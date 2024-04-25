#!/usr/bin/env python3

import argparse
import sys

import pandas as pd

from extract_parquet_results import QtlLocusVariantFilter, QtlResultProcessor, QtlGeneFilter


def calculate_ld(input_file, output_file, variant_filters, gene_filter):
    result_processor = QtlResultProcessor(
        input_file, gene_filter=gene_filter)

    result_processor.variant_filters = variant_filters

    df = result_processor.extract(
        add={'z_score'}) 

    # pivot the DataFrame to get a matrix of z-scores
    matrix = df.pivot(index='phenotype', columns='variant', values='z_score')

    print(f"Running matrix.corr() on a {matrix.shape[0]} by {matrix.shape[1]} matrix... (this might take a while)")
    # calculate the pairwise correlations between variants
    corr_matrix = matrix.corr()

    # write the list of uncorrelated genes to a file
    corr_matrix.to_csv(output_file)
    print("Done!")
    print("Closing output file '{}'".format(output_file))


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description="Analyze HASE output parquet files.")

    parser.add_argument('-i', '--input-file', required=True,
                        help="""One or multiple input parquet arrays. 
                                  Usage of wildcards is supported but then the argument has to be quoted.""")
    parser.add_argument('-o', '--output-prefix', type=str,
                        required=True,
                        help="Path to tab-separated output file.")
    parser.add_argument('-B', '--bed-file', required=True, default=None,
                        help="""Bed file with loci to extract""")
    parser.add_argument('-G', '--genes-file', required = True, default = None,
                        help = """File with the list of phenotypes to include.""")
    parser.add_argument('-r', '--variant-reference', required=True,
                        help="Reference for variants. Has to be gzipped and space-delimited.")

    args = parser.parse_args(argv[1:])

    print(args)

    variant_reference = (
        pd.read_csv(args.variant_reference, sep=' ')
        .drop(["allele1", "allele2"], axis=1)
        .rename({"ID": "variant", "bp": "bp", "CHR": "chromosome", "str_allele1": "a1", "str_allele2": "a2"},
                axis=1))
    print(variant_reference.head())

    print("Using loci file '%s' to filter on variants." % args.variant_reference)
    loci = pd.read_csv(args.bed_file, sep="\t", header=None, names=["chromosome", "start", "stop"])

    print("Using variants file '%s' to filter on variants." % args.genes_file)
    qtl_gene_filter = QtlGeneFilter.from_path(args.genes_file)

    for i, (index, row) in enumerate(loci.iterrows()):
        print("Starting export for locus {}/{}".format(i + 1, loci.shape[0]))

        chromosome = row["chromosome"]
        start = row["start"]
        stop = row["stop"]

        locus_filter = QtlLocusVariantFilter.from_locus(
            chromosome, start, stop, variant_reference)

        calculate_ld(args.input_file,
                     "{}.{}_{}_{}.csv.gz".format(args.output_prefix, chromosome, start, stop),
                     [locus_filter], qtl_gene_filter)

    return 0


if __name__ == '__main__':
    sys.exit(main())
