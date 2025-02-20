#!/usr/bin/env python3

"""
Created:      29/11/2024
Author:       C.A. (Robert) Warmerdam

Copyright (C) 2024 C.A. Warmerdam

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

A copy of the GNU General Public License can be found in the LICENSE file in the
root directory of this source tree. If not, see <https://www.gnu.org/licenses/>.
"""

import sys
import argparse
import pandas as pd
import pyarrow.parquet as pq
import pyarrow.feather as ft
import numpy as np
import time
import ld_utils
from carma_gentropy import CARMA


def extract_summary_stats(empirical_dataset, gene_cluster, variant_reference):
    gene = gene_cluster["gene"].unique()[0]
    locus_chromosome = gene_cluster["chromosome"].unique()[0]
    gene_summary_stats_list = list()
    print(f"Extracting summary statistics for gene {gene}, {gene_cluster.shape[0]} loci")
    for start, end in zip(gene_cluster["start"], gene_cluster["end"]):
        print(f"{locus_chromosome}:{start}-{end}")
        # Get variants in gene region
        gene_locus_variant_reference = variant_reference[
            (variant_reference["chromosome"] == locus_chromosome) &
            (variant_reference["bp"].between(start, end))
            ]

        gene_variant_index_start = gene_locus_variant_reference["variant_index"].min()
        gene_variant_index_end = gene_locus_variant_reference["variant_index"].max()

        gene_summary_stats_list.append(pq.ParquetDataset(
            empirical_dataset,
            filters=[[
                ("phenotype", "==", gene),
                ("variant_index", ">=", gene_variant_index_start),
                ("variant_index", "<=", gene_variant_index_end)]]).read().to_pandas())
        print(f"Found {gene_summary_stats_list[-1].shape[0]} variants")
    gene_summary_stats = pd.concat(gene_summary_stats_list, axis=0)
    print(f"Total number of variants for {gene}: {gene_summary_stats.shape[0]} variants")
    return gene_summary_stats


def filter_locus(
        empirical_dataset, ld_calculator, locus_bed, variant_reference,
        min_sample_size_prop=0.8, max_i_squared=40, min_i_squared_to_inspect=40, normalize_sumstats=True,
        debug=False, dry_run=False
):

    ld_matrix, variant_order = ld_calculator.calculate_ld_non_continuous(
        locus_bed)
    variant_index_start = min(variant_order)
    variant_index_end = max(variant_order)
    variant_indices = dict(zip(variant_order, range(len(variant_order))))

    if debug == 'locus-wise':
        ld_matrix.to_csv(f"LD_{variant_index_start}-{variant_index_end}.csv")
    elif debug == 'write-susie':
        ft.write_feather(
            pd.DataFrame(ld_matrix, index=variant_order, columns=variant_order),
            f"LD_{variant_index_start}-{variant_index_end}.feather"
        )

    # Fine-mapping the locus
    fine_mapping_results = list()

    for index, gene_cluster in locus_bed.group_by(["gene", "gene_cluster"]):
        gene_summary_stats = extract_summary_stats(empirical_dataset, gene_cluster, variant_reference)

        # Apply sample size filtering
        sample_size_threshold = gene_summary_stats["sample_size"].max() * min_sample_size_prop
        gene_summary_stats = gene_summary_stats[
            (gene_summary_stats["sample_size"] >= sample_size_threshold) &
            (gene_summary_stats["i_squared"] <= max_i_squared)
            ]

        gene_summary_stats["Z"] = gene_summary_stats["beta"] / gene_summary_stats["standard_error"]
        # Normalize summary statistics
        if normalize_sumstats:
            gene_summary_stats["beta"] = gene_summary_stats["Z"] / np.sqrt(
                gene_summary_stats["sample_size"] + gene_summary_stats["Z"] ** 2)
            gene_summary_stats["standard_error"] = 1 / np.sqrt(
                gene_summary_stats["sample_size"] + gene_summary_stats["Z"] ** 2)

        # Filter variant order
        variant_order_filtered = [v for v in variant_order if v in gene_summary_stats["variant_index"].values]
        gene_summary_stats = gene_summary_stats.set_index("variant_index").loc[variant_order_filtered].reset_index()
        gene_summary_stats['outliers'] = False

        print(gene_summary_stats)

        # Perform fine-mapping if not a dry run
        if len(gene_summary_stats) > 0 and all(
                gene_summary_stats["variant_index"] == variant_order_filtered) and not dry_run:
            # Fine-mapping step (to be implemented)
            ld_variant_indices = np.array([variant_indices[variant_name] for variant_name in variant_order_filtered])
            print(f"Starting CARMA outlier detection with {len(ld_variant_indices)} variants!")
            carma_out = CARMA.CARMA_spike_slab_noEM(
                gene_summary_stats['Z'], ld_matrix[np.ix_(ld_variant_indices, ld_variant_indices)], all_iter=1)
            outliers = carma_out["Outliers"]
            print(f"Found {len(outliers)} outliers!")

            gene_summary_stats.loc[np.array(outliers), ['outliers']] = True
            fine_mapping_results.append(gene_summary_stats)

    return pd.concat(fine_mapping_results, ignore_index=True) if fine_mapping_results else None


def main(argv=None):
    if argv is None:
        argv = sys.argv

        # Argument Parser
    parser = argparse.ArgumentParser(description='Run fine-mapping over loci.')
    parser.add_argument("--ld", required=True, help="Path to the directory containing ld data.")
    parser.add_argument("--ld-type", required=True, choices=["gene-set", "pcs", "dosages"], help="Type of LD")
    parser.add_argument("--empirical", required=True, help="Path to the directory containing empirical data.")
    parser.add_argument("--variant-reference", required=True, help="Path to the variant reference file.")
    parser.add_argument("--uncorrelated-genes", required=False, help="File containing uncorrelated genes.")
    parser.add_argument("--max-i2", required=False, default=40, type=float, help="Maximum i2")
    parser.add_argument("--inspect-min-i2", required=False, default=40, type=float, help="Minimum i2 value to inspect for outliers")
    parser.add_argument("--min-n-prop", required=False, default=0.8, type=float,
                        help="Minimum sample size proportion compared to maximum in locus")
    parser.add_argument("--no-adjust-stats", required=False, action='store_true',
                        help="Flag to disable adjusting betas and standard errors for sample size")
    parser.add_argument("--bed-files", required=True, nargs='+', help="Space-separated list of BED files.")
    parser.add_argument("--debug", required=False, default='off', help="Enables debugging mode.",
                        choices=['off', 'locus-wise', 'write-susie'])
    parser.add_argument("--dry-run", required=False, action='store_true',
                        help="Runs code without actually running fine-mapping")

    args = parser.parse_args()

    # Print argument values
    print(f"LD data directory: {args.ld}")
    print(f"LD type: {args.ld_type}")
    print(f"Empirical data directory: {args.empirical}")
    print(f"Variant reference file: {args.variant_reference}")
    print(f"Uncorrelated genes file: {args.uncorrelated_genes}")
    print(f"BED files: {', '.join(args.bed_files)}")
    print(f"Max I2: {args.max_i2}")
    print(f"Min I2 to filter: {args.inspect_min_i2}")
    print(f"Min sample size prop: {args.min_n_prop}")
    print(f"No adjust stats: {args.no_adjust_stats}")
    print(f"Debugging: {args.debug}")
    print(f"Dry-run: {args.dry_run}")

    # Load input files
    variant_reference = pq.read_table(args.variant_reference).to_pandas()
    uncorrelated_genes = pd.read_csv(args.uncorrelated_genes, header=None)[0].tolist()

    normalize_sumstats = not args.no_adjust_stats
    min_sample_size_prop = args.min_n_prop
    max_i_squared = args.max_i2
    min_i_squared_to_inspect = args.inspect_min_i2
    dry_run = args.dry_run
    debug = args.debug

    # Load empirical dataset
    empirical_dataset = args.empirical

    # LD Function Definition
    if args.ld_type == "gene-set":
        ld_calculator = ld_utils.LdCalculator(
            parquet_dataset=args.ld, observations=uncorrelated_genes, variant_reference=variant_reference)
    elif args.ld_type == "pcs":
        ld_calculator = ld_utils.LdCalculatorWide(
            parquet_dataset=args.ld, observations=None, variant_reference=variant_reference
        )
    else:
        raise NotImplementedError("Cannot find implementation for other LD panel type.")

    # Fine-mapping process
    filtering_results_per_locus = []

    for bed_file in args.bed_files:
        locus_bed = pd.read_csv(bed_file, sep="\t", names=["chromosome", "start", "end", "gene", "gene_cluster", "cluster"])
        print(
            f"Starting finemapping in {locus_bed['chromosome'].iloc[0]}:{locus_bed['start'].min()}-{locus_bed['end'].max()} ({len(locus_bed)} genes)")

        filter_output = filter_locus(
            empirical_dataset=empirical_dataset,
            ld_calculator=ld_calculator,
            locus_bed=locus_bed,
            variant_reference=variant_reference,
            min_sample_size_prop=min_sample_size_prop,
            max_i_squared=max_i_squared,
            min_i_squared_to_inspect=min_i_squared_to_inspect,
            normalize_sumstats=normalize_sumstats,
            debug=debug,
            dry_run=dry_run
        )

        filtering_results_per_locus.append(filter_output)

    # Combine results
    combined_results = pd.concat(
        filtering_results_per_locus,
        ignore_index=True) if filtering_results_per_locus else None

    # Save results
    if not dry_run and combined_results is not None:
        combined_results.to_csv("carma.results.tsv", sep="\t", index=False)
    if debug == 'write-susie':
        ft.write_feather(
            combined_results[~combined_results.outliers],
            f"susie_input.feather"
        )


if __name__ == "__main__":
    sys.exit(main())
