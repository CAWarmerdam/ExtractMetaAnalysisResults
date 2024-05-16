#!/usr/bin/env python3

import argparse
import sys

import numpy as np
import pandas as pd


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(description = "Analyze HASE output parquet files.")

    parser.add_argument('-P', '--gene-variant-pairs-file', required=False, default=None, help="""Gene variant pairs to extract""")
    parser.add_argument('-k', '--k-genes', type=int)

    args = parser.parse_args(argv[1:])

    pairs = pd.read_csv(args.gene_variant_pairs_file, delimiter='\t', header=1, names=["variant", "gene"])
    print(pairs)
    pairs_grouped = pairs.groupby('gene')
    chunk = 0
    n_genes = 0

    for gene, group in pairs_grouped:
        n_genes += 1
        output_file = f'output_chunk_{chunk}.csv'
        print(output_file)
        group.to_csv(output_file, sep="\t", header=False, index=False, mode='a')
        if n_genes == args.k_genes:
            n_genes = 0
            chunk += 1

    return 0


if __name__ == '__main__':
    sys.exit(main())
