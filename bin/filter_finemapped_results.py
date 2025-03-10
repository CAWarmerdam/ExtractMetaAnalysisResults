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
import pandas as pd


def process_and_save_dataframe(df, output_folder):
    # Ensure required columns exist
    if 'SusieRss_lambda' not in df.columns:
        df['SusieRss_lambda'] = None
    if 'SusieRss_pip' not in df.columns:
        df['SusieRss_pip'] = None

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


def main(args):
    df_list_pass = list()
    for path in args.inputPath:
        df = pd.read_csv(path, sep="\t")
        df_unfiltered = df.copy()
        if 'naive' in args.strategy:
            df = df[df['SusieRss_pip'] > 0.9]
        if 'lbf' in args.strategy:
            lbf_cols = [col for col in df.columns if col.startswith('lbf_cs_')]
            df["max_lbf"] = df[lbf_cols].max(axis=1)
            df = df[(df["max_lbf"] > 2) & df["SusieRss_CS"].notna()]
        if 'cs' in args.strategy:
            df = df[df["SusieRss_CS"].notna()]
        df_list_pass.append(df)

        if 'write_failed_loci' in args.strategy:
            process_and_save_dataframe(df_unfiltered, args.failedOutputPath)

    df_concat_pass = pd.concat(df_list_pass, sort=True)
    df_concat_pass.to_csv(args.outputPath, sep="\t", index=False)
    return


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--inputPath", type=str, required=True, help="Help goes here", nargs='+')
    parser.add_argument("-s", "--strategy", type=str, required=True, help="Help goes here", nargs='+')
    parser.add_argument("-o", "--outputPath", type=str, required=True, help="Help goes here")
    parser.add_argument("-f", "--failedOutputPath", type=str, required=False, help="Path to output parquet dataset")
    args = parser.parse_args()

    main(args)
