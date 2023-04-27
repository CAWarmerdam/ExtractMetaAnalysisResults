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
    parser.add_argument('-g', '--genes', required=False, default=None, nargs = '+',
                        help = """Individual phenotype IDs specified and separated by space.""")
    parser.add_argument('-r', '--variant-reference', dest='variant_reference', required=True,
                        help='Path to the table containing all SNPs from a reference panel')
    parser.add_argument('-l', '--loci', required = True, default = None,
                        help = """A bed file that contains one or more loci for which to calculate LD""")

    args = parser.parse_args(argv[1:])

    variant_reference = (
        pd.read_csv(args.variant_reference, sep = ' ')
        .drop(["allele1", "allele2"], axis=1)
        .rename({"ID": "variant", "bp": "bp", "CHR": "chromosome", "str_allele1": "a1", "str_allele2": "a2"}, axis=1))

    qtl_gene_filter = None

    if args.genes is not None:
        print("Provided %d genes for filtering." % len(args.genes))
        qtl_gene_filter = QtlGeneFilter.from_list(args.genes)

    loci = pd.read_csv(args.loci, sep=" ", header=None, names=["chromosome", "start", "stop", "name"])
    for index, row in loci.iterrows():

        locus = row["name"]
        chromosome = row["chromosome"]
        start = row["start"]
        stop = row["stop"]

        locus_filter = QtlLocusVariantFilter.from_locus(
            chromosome, start, stop, variant_reference)

        result_processor = QtlResultProcessor(
            args.input_file, gene_filter=qtl_gene_filter)
        result_processor.variant_filters = [locus_filter]

        # Get dataframe with z-scores
        df = result_processor.extract(cols="z_score")

        # pivot the DataFrame to get a matrix of z-scores
        matrix = df.pivot(index='variant', columns='phenotype', values='z_score')

        # calculate the pairwise correlations between variants
        corr_matrix = matrix.corr()

        # write the list of uncorrelated genes to a file
        corr_matrix.to_csv("{prefix}.{chrom}_{start}-{stop}_{name}".format(
            prefix=args.output_prefix, chrom=chromosome, start=start, stop=stop, name=locus))

    return 0


if __name__ == '__main__':
    sys.exit(main())
