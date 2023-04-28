#!/usr/bin/env python3

import argparse
import sys

import pandas as pd
import pyarrow as pa

from extract_parquet_results import QtlResultProcessor, QtlGeneFilter, \
    QtlLocusVariantFilter

SCHEMA = pa.schema([("variant", pa.string()), ("beta", pa.float64()),
                    ("standard_error", pa.float64()), ("i_squared", pa.float64()),
                    ("sample_size", pa.float64())])


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description = "Extract HASE output for specified loci, and calculate LD matrix "
                                                   "for these loci.")

    parser.add_argument('-i', '--input-file', required = True,
                        help = """Input parquet dataset""")
    parser.add_argument('-o', '--output-prefix', type = str,
                        required = True,
                        help = "Output prefix which to use for writing LD data.")
    parser.add_argument('-G', '--genes-file', required = False, default = None,
                        help = """File with the list of phenotypes to include.""")
    parser.add_argument('-r', '--variant-reference', dest='variant_reference', required=True,
                        help='Path to the table containing all SNPs from a reference panel')
    parser.add_argument('-l', '--loci', required = True, default = None,
                        help = """A bed file that contains one or more loci for which to calculate LD""")

    args = parser.parse_args(argv[1:])
    print(args)

    variant_reference = (
        pd.read_csv(args.variant_reference, sep = ' ')
        .drop(["allele1", "allele2"], axis=1)
        .rename({"ID": "variant", "bp": "bp", "CHR": "chromosome", "str_allele1": "a1", "str_allele2": "a2"}, axis=1))

    qtl_gene_filter = None

    if args.genes_file is not None:
        print("Using variants file '%s' to filter on variants." % args.genes_file)
        qtl_gene_filter = QtlGeneFilter.from_path(args.genes_file)

    loci = pd.read_csv(args.loci, sep="\t", header=None, names=["chromosome", "start", "stop", "name"])
    print(loci)

    for index, row in loci.iterrows():
        print(row)

        locus = row["name"].split(",")
        print(locus)
        chromosome = row["chromosome"]
        start = row["start"]
        stop = row["stop"]

        locus_filter = QtlLocusVariantFilter.from_locus(
            chromosome, start, stop, variant_reference)

        result_processor = QtlResultProcessor(
            args.input_file, gene_filter=qtl_gene_filter)
        result_processor.variant_filters = [locus_filter]

        # Get dataframe with z-scores
        df = result_processor.extract(cols={"z_score",})

        # pivot the DataFrame to get a matrix of z-scores
        matrix = df.pivot(index='phenotype', columns='variant', values='z_score')

        # calculate the pairwise correlations between variants
        corr_matrix = matrix.corr()

        # write the list of uncorrelated genes to a file
        corr_matrix.to_csv("{prefix}.{chrom}_{start}-{stop}.csv".format(
            prefix=args.output_prefix, chrom=chromosome, start=start, stop=stop))

    return 0


if __name__ == '__main__':
    sys.exit(main())
