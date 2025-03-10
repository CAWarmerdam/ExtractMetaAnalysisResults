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
import RSparsePro
from run_carma_light_over_loci import extract_summary_stats


def finemap_locus(
        empirical_dataset, ld_calculator, locus_bed, variant_reference,
        min_sample_size_prop=0.8, max_i_squared=40, normalize_sumstats=True,
        debug=False, dry_run=False, nCS=20, ft=None):

    ld_matrix, variant_order, variant_index_start, variant_index_end = ld_calculator.calculate_ld_non_continuous(
        locus_bed)
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

    for index, gene_cluster in locus_bed.groupby(["gene", "gene_cluster"]):
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
        if len(gene_summary_stats) > 0 and all(gene_summary_stats["variant_index"] == variant_order_filtered) and not dry_run:
            # Fine-mapping step (to be implemented)
            maxite = 100
            eps = 1e5
            ubound = 100000
            cthres = 0.95
            eincre = 1.5
            minldthres = 0.7
            maxldthres = 0.2
            varemax = 100.0
            varemin = 1e-3
            ld_variant_indices = np.array([variant_indices[variant_name] for variant_name in variant_order_filtered])
            eff, eff_gamma, eff_mu, PIP, ztilde = RSparsePro.adaptive_train(
                gene_summary_stats["Z"],
                ld_matrix[np.ix_(ld_variant_indices, ld_variant_indices)],
                nCS, maxite, eps, ubound, cthres, minldthres, maxldthres, eincre, varemax, varemin)
            gene_summary_stats['pip'] = PIP
            gene_summary_stats['z_estimated'] = ztilde
            gene_summary_stats['cs'] = 0
            print(f'The maximum Z score in the window: {gene_summary_stats["Z"].max()}')
            for e in eff:
                mcs_idx = [gene_summary_stats['variant_index'][j] for j in eff[e]]
                print(f'The {e}-th effect group contains effective variants:')
                print(f'causal variants: {mcs_idx}')
                print(f'variant probabilities for this effect group: {eff_gamma[e]}')
                print(f'zscore for this effect group: {eff_mu[e]}\n')
                gene_summary_stats.iloc[eff[e], gene_summary_stats.columns.get_loc('cs')] = e+1
        else:
            gene_summary_stats['pip'] = np.nan
            gene_summary_stats['z_estimated'] = np.nan
            gene_summary_stats['cs'] = np.nan

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
    parser.add_argument("--min-n-prop", required=False, default=0.8, type=float, help="Minimum sample size proportion compared to maximum in locus")
    parser.add_argument("--no-adjust-stats", required=False, action='store_true', help="Flag to disable adjusting betas and standard errors for sample size")
    parser.add_argument("--bed-files", required=True, nargs='+', help="Space-separated list of BED files.")
    parser.add_argument("--debug", required=False, default='off', help="Enables debugging mode.", choices=['off', 'locus-wise'])
    parser.add_argument("--dry-run", required=False, action='store_true', help="Runs code without actually running fine-mapping")

    args = parser.parse_args()

    # Print argument values
    print(f"LD data directory: {args.ld}")
    print(f"LD type: {args.ld_type}")
    print(f"Empirical data directory: {args.empirical}")
    print(f"Variant reference file: {args.variant_reference}")
    print(f"Uncorrelated genes file: {args.uncorrelated_genes}")
    print(f"BED files: {', '.join(args.bed_files)}")
    print(f"Max I2: {args.max_i2}")
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
    fine_mapping_results_per_locus = []

    for bed_file in args.bed_files:
        locus_bed = pd.read_csv(bed_file, sep="\t", names=["chromosome", "start", "end", "gene", "gene_cluster", "cluster"])
        print(f"Starting finemapping in {locus_bed['chromosome'].iloc[0]}:{locus_bed['start'].min()}-{locus_bed['end'].max()} ({len(locus_bed)} genes)")

        fine_mapping_output = finemap_locus(
            empirical_dataset=empirical_dataset,
            ld_calculator=ld_calculator,
            locus_bed=locus_bed,
            variant_reference=variant_reference,
            min_sample_size_prop=min_sample_size_prop,
            max_i_squared=max_i_squared,
            normalize_sumstats=normalize_sumstats,
            debug=debug,
            dry_run=dry_run
        )

        fine_mapping_results_per_locus.append(fine_mapping_output)

    # Combine results
    combined_results = pd.concat(fine_mapping_results_per_locus, ignore_index=True) if fine_mapping_results_per_locus else None

    # Save results
    if not dry_run and combined_results is not None:
        combined_results.to_csv("finemapped.results.tsv", sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())
