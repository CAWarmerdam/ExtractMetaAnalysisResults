#!/usr/bin/env python2

import argparse
import re
import sys
from abc import ABCMeta, abstractmethod

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from scipy import stats

from bin.extract_parquet_results import QtlResultProcessor, QtlGeneFilter, QtlVariantFilter, QtlPThresholdFilter

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
    parser.add_argument('-r', '--variant-reference', dest='variant_reference', required=True,
                        help='Path to the table containing all SNPs from a reference panel')
    parser.add_argument('-l', '--loci', required = True, default = None,
                        help = """A bed file that contains one or more loci for which to calculate LD""")

    args = parser.parse_args(argv)

    result_processor = QtlResultProcessor(
        args.input_file)

    # Get dataframe with z-scores
    df = result_processor.extract(cols="z-score")

    # pivot the DataFrame to get a matrix of z-scores
    matrix = df.pivot(index='variant', columns='phenotype', values='z_score')

    # calculate the pairwise correlations between genes
    corr_matrix = matrix.corr()



    return 0
