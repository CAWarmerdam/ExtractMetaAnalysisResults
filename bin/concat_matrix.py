#!/usr/bin/env python3


import sys
import argparse
import pandas as pd
import numpy as np
import numpy.ma as ma


def main(argv=None):
    if argv is None:
        argv = sys.argv

    # parse command-line arguments
    parser = argparse.ArgumentParser(description='Calculate correlations between z-scores and find uncorrelated genes')
    parser.add_argument('--input', nargs='+')
    parser.add_argument('--output')
    args = parser.parse_args()

    concatenated = pd.concat([pd.read_csv(input_matrix, sep="\t", index_col=0, header=0) for input_matrix in args.input], axis=1)
    print(concatenated.head())
    concatenated.to_csv(args.output, sep="\t", index_label="variant")
 
    return 0


if __name__ == '__main__':
    sys.exit(main())
