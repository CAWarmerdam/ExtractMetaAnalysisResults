#!/usr/bin/env python3

import argparse
import sys

import pandas as pd
import pyarrow as pa


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description = "Analyze HASE output parquet files.")

    parser.add_argument('-i', '--input-file', required = True,
                        help = """One or multiple input parquet arrays. 
                                  Usage of wildcards is supported but then the argument has to be quoted.""")
    parser.add_argument('-o', '--output-prefix', type = str,
                        required = True,
                        help = "Path to tab-separated output file.")
    parser.add_argument('-B', '--bed-file', required=False, default=None,
                        help = """Bed file with loci to extract""")
    parser.add_argument('-r', '--variant-reference', required = False,
                        help = "Reference for variants. Has to be gzipped and space-delimited.")

    args = parser.parse_args(argv[1:])

    print(args)

    variant_reference = None
    variants_list = None
    variant_filters = None
    loci = None

    if args.variant_reference is not None:
        variant_reference = (
            pd.read_csv(args.variant_reference, sep = ' ')
            .drop(["allele1", "allele2"], axis=1)
            .rename({"ID": "variant", "bp": "bp", "CHR": "chromosome", "str_allele1": "a1", "str_allele2": "a2"}, axis=1))
        print(variant_reference.head())
    if args.bed_file is not None:
        print("Using loci file '%s' to filter on variants." % args.variant_reference)
        if variants_list is not None:
            print("Variant filter already defined. Skipping...")
        loci = pd.read_csv(args.bed_file, sep="\t", header=None, names=["chromosome", "start", "stop", "name"])

    if variants_list is not None:
        if variant_reference is None:
            parser.error("Cannot subset on variants without variant reference")
        variant_selection = variant_reference.loc[variant_reference.loc[:, "variant"].isin(variants_list), :]
        variant_dictionary = variant_selection.groupby('chromosome')['variant'].apply(list).to_dict()
        variant_filters = [QtlLocusVariantFilter(chromosome, variants) for chromosome, variants in variant_dictionary.items()]

    else:
        for i, (index, row) in enumerate(loci.iterrows()):
            print("Starting export for locus {}/{}".format(i+1, loci.shape[0]))

            locus = row["name"].split(",")
            chromosome = row["chromosome"]
            start = row["start"]
            stop = row["stop"]

            locus_filter = QtlLocusVariantFilter.from_locus(
                chromosome, start, stop, variant_reference)

            if qtl_gene_filter is None:
                qtl_gene_filter = QtlGeneFilter.from_list(locus)

            output_file = "{}.{}_{}-{}.out.csv".format(args.output_prefix, chromosome, start, stop)

            export_write(args.input_file, output_file,
                         qtl_gene_filter, [locus_filter],
                         args.column_specifications, args.p_thresh)

    for i, (index, row) in enumerate(loci.iterrows()):
        print("Starting export for locus {}/{}".format(i+1, loci.shape[0]))

        chromosome = row["chromosome"]
        start = row["start"]
        stop = row["stop"]

        locus_file = "{}.{}_{}-{}.out.csv".format(args.output_prefix, chromosome, start, stop)

        df = pd.read_csv(locus_file)

        # pivot the DataFrame to get a matrix of z-scores
        matrix = df.pivot(index='phenotype', columns='variant', values='z_score')

        # calculate the pairwise correlations between variants
        corr_matrix = matrix.corr()

        # write the list of uncorrelated genes to a file
        corr_matrix.to_csv("{prefix}.{chrom}_{start}-{stop}.csv.gz".format(
            prefix=args.output_prefix, chrom=chromosome, start=start, stop=stop))

    return 0


if __name__ == '__main__':
    sys.exit(main())
