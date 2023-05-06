#!/usr/bin/env python3

import argparse
import sys

import pandas as pd
import pyarrow as pa


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description = "Extract HASE output for specified loci, and calculate LD matrix "
                                                   "for these loci.")

    parser.add_argument('-i', '--input-prefix', required = True,
                        help = """Input parquet dataset""")
    parser.add_argument('-o', '--output-prefix', type = str,
                        required = True,
                        help = "Output prefix which to use for writing LD data.")
    parser.add_argument('-l', '--loci', required = False, default = None,
                        help = """File with the list of phenotypes to include.""")

    args = parser.parse_args(argv[1:])
    print(args)

    loci = pd.read_csv(args.loci, sep="\t", header=None, names=["chromosome", "start", "stop", "name"])
    print(loci)

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
