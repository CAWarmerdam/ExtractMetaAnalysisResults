#!/usr/bin/env python3
"""
<A single line describing this program goes here.>

MIT License

Copyright (c) 2022 Tijs van Lieshout

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

Uses:
<The terminal interactions with this script go here>
"""

# Metadata
__title__ = "Template for a CLI python script"
__author__ = "Tijs van Lieshout"
__created__ = "2022-05-04"
__updated__ = "2022-05-27"
__maintainer__ = "Tijs van Lieshout"
__email__ = "t.van.lieshout@umcg.nl"
__version__ = 0.2
__license__ = "GPLv3"
__description__ = f"""{__title__} is a python script created on {__created__} by {__author__}.
                      Last update (version {__version__}) was on {__updated__} by {__maintainer__}.
                      Under license {__license__} please contact {__email__} for any questions."""

# Imports
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm
from tqdm import tqdm


def process_and_save_failed_genes(df, output_folder):
    # Group by 'phenotype' and 'SusieRss_lambda'
    grouped = df.groupby(['phenotype', 'SusieRss_lambda'])

    # Filter groups where all 'SusieRss_pip' values are NaN
    filtered_groups = {name: group for name, group in grouped if group['SusieRss_pip'].isna().all()}

    # Keep only selected columns
    selected_columns = ['variant_index', 'standard_error', 'sample_size', 'phenotype', 'beta', 'i_squared']

    # Ensure output folder exists
    os.makedirs(output_folder, exist_ok=True)

    # Write each group to a partitioned Parquet file
    for (phenotype, _), group in filtered_groups.items():
        min_variant_index = group['variant_index'].min()
        max_variant_index = group['variant_index'].max()
        output_path = os.path.join(
            output_folder, f"phenotype={phenotype}_loc={min_variant_index}-{max_variant_index}.feather")

        group[selected_columns].to_feather(output_path)
        print(f"Saved: {output_path}")


def summarize_loci(df):
    summaries = []
    grouped = df.groupby(['phenotype', 'SusieRss_lambda'], dropna=False)
    for (phenotype, lambda_value), group in grouped:
        lbf_cols = [col for col in df.columns if col.startswith('lbf_cs_')]
        if lbf_cols:
            group["max_lbf"] = group[lbf_cols].max(axis=1)

        chi_squared_values = np.power(group['beta'] / group['standard_error'], 2)

        group['p_value'] = 2 * (1 - norm.cdf(np.abs(group['beta'] / group['standard_error'])))

        summary = {
            "phenotype": phenotype,
            "SusieRss_lambda": lambda_value,
            "num_snps": len(group),
            "converged": group['SusieRss_pip'].notna().all(),
            "min_variant_index": group['variant_index'].min(),
            "max_variant_index": group['variant_index'].max(),
            "num_unique_SusieRss_CS": group['SusieRss_CS'].nunique(dropna=True),
            "num_SusieRss_CS_lbf": group[(group['max_lbf'] > 2) & group['SusieRss_CS'].notna() & group['p_value'] < 1e-5].shape[0],
            "min_sample_size": group['sample_size'].min(),
            "max_sample_size": group['sample_size'].max(),
            "max_chi2": chi_squared_values.max(),
            "mean_chi2": chi_squared_values.mean(),
            "max_i2": group['i_squared'].max(),
            "median_i2": group['i_squared'].median(),
            "mean_i2": group['i_squared'].mean()
        }
        summaries.append(summary)
    return pd.DataFrame(summaries)


def parse_finemapping_output(args):
    df_list_pass = list()
    df_list_summary = list()
    for path in tqdm(args.inputPath):
        df = pd.read_csv(path, sep="\t")
        # Ensure required columns exist
        if 'SusieRss_lambda' not in df.columns:
            df['SusieRss_lambda'] = None
        if 'SusieRss_pip' not in df.columns:
            df['SusieRss_pip'] = None
        if 'SusieRss_CS' not in df.columns:
            df['SusieRss_CS'] = None
        df_unfiltered = df.copy()
        if 'none' not in args.strategy:
            if 'naive' in args.strategy:
                df = df[df['SusieRss_pip'] > 0.9]
            if 'lbf' in args.strategy:
                lbf_cols = [col for col in df.columns if col.startswith('lbf_cs_')]
                df["max_lbf"] = df[lbf_cols].max(axis=1)
                df = df[(df["max_lbf"] > 2) & df["SusieRss_CS"].notna()]
            if 'cs' in args.strategy:
                df = df[df["SusieRss_CS"].notna()]
            df_list_pass.append(df)

        if args.failedOutputPath is not None:
            process_and_save_failed_genes(df_unfiltered, args.failedOutputPath)
        if args.summaryOutputPath is not None:
            df_list_summary.append(summarize_loci(df_unfiltered))

    if len(df_list_pass) > 0:
        df_concat_pass = pd.concat(df_list_pass, sort=True)
        df_concat_pass.to_csv(args.outputPath, sep="\t", index=False)

    if len(df_list_summary) > 0:
        df_concat_summary = pd.concat(df_list_summary, sort=True)
        df_concat_summary.to_csv(args.summaryOutputPath, sep="\t", index=False)


def main(argv=None):
    if argv is None:
        argv = sys.argv

    parser = argparse.ArgumentParser(argv[1:])

    parser.add_argument("-i", "--inputPath",
                        type=str, required=True, help="Help goes here", nargs='+')
    parser.add_argument("-s", "--strategy",
                        type=str, required=True, help="Help goes here", nargs='+')
    parser.add_argument("-o", "--outputPath",
                        type=str, required=True, help="Help goes here")
    parser.add_argument("-f", "--failedOutputPath",
                        type=str, required=False, help="Path to output parquet dataset", default=None)
    parser.add_argument("-S", "--summaryOutputPath",
                        type=str, required=False, help="Path to write summary to", default=None)
    args = parser.parse_args()

    parse_finemapping_output(args)
    return


if __name__ == '__main__':
    sys.exit(main())
