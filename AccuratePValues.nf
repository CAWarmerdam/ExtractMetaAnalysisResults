#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { ExtractSignificantResults; AnnotateLoci; IntersectLoci } from './modules/CollectSignificantLoci'
include { CalculateAccuratePermutationPValues; GetUncorrelatedVariantsWithLdDataset } from './modules/CalculateAccuratePermutationPValues'
include { CalculateLdMatrix; UncorrelatedGenes } from './modules/CalculateLdMatrix'
include { GetUncorrelatedVariants } from './modules/UncorrelatedVariants'
include { CalculateZScores } from './modules/CalculateZScores'


def helpmessage() {

log.info"""

HASE output analyzer v${workflow.manifest.version}
==============================================
Pipeline for parallelized extraction and filtering of the raw HASE results.

This pipeline is used to extract subsets of results from the HASE results (numerous large .parquet files).

Usage:

nextflow run ExtractHaseResults.nf \
--empirical '/inputfolder/' \
--permuted '/outputfile/' \
--genes '/phenotypes.txt' \
--ld-dataset '/dataset/' \
--output '/output/'


Mandatory arguments:
--empirical           Path to the folder with HASE result .parquet files.
--permuted            Path to where the database should be written
--genes               Path to a file with all unique genes
--maf-table           Path to table with maf per variant
--ld-dataset          Path to LD dataset
--output         	  Path to outputfolder

""".stripIndent()

}

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.empirical).collect().set { empirical_parquet_ch }
Channel.fromPath(params.permuted).collect().set { permuted_parquet_ch }
Channel.fromPath(params.reference_data).set { reference_bcf_files_ch }
Channel.fromPath(params.genes).splitText().set { genes_ch }
Channel.fromPath(params.genome_reference).collect().set { genome_ref_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.gene_reference).collect().set { gene_reference_ch }

log.info """=======================================================
HASE output analyzer v${workflow.manifest.version}"
======================================================="""
def summary = [:]
summary['Pipeline Version']                         = workflow.manifest.version
summary['Current user']                             = "$USER"
summary['Current home']                             = "$HOME"
summary['Current path']                             = "$PWD"
summary['Working dir']                              = workflow.workDir
summary['Script dir']                               = workflow.projectDir
summary['Config Profile']                           = workflow.profile
summary['Container Engine']                         = workflow.containerEngine
if(workflow.containerEngine) summary['Container']   = workflow.container
summary['Empirical eQTLs']                          = params.empirical
summary['Permuted eQTLs']                           = params.permuted
summary['Reference data']                           = params.reference_data
summary['Genome reference']                         = params.genome_reference
summary['Variant reference']                        = params.variant_reference
summary['Gene reference']                           = params.gene_reference
summary['Gene list']                                = params.genes

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Procedure:
// Generate breakpoints to get equal numbers of variants per bin
// Get uncorrelated variants per maf bin
// Using the uncorrelated variants, do accurate permutation p-value calculation

// workflow ACCURATE_P_VALUES {
//     take:
//         mafTable
//         empiricalResults
//         permutedResults
//         breakPoints
//
//     main:
//         GetBreakpoints(mafTable.collect())
//         uncorrelatedVariants = GetUncorrelatedVariants(GetBreakpoints.splitText())
//         CalculateAccuratePermutationPValues(
//             empiricalResults, permutedResults,
//             mafTable, breakpoints, uncorrelatedVariants, breakPoints, uncorrelatedVariants, andersonDarlingTable)
//
//     emit:
//         CalculateAccuratePermutationPValues.out
// }

workflow CALCULATE_LD {
        // Obtain a list of uncorrelated variants
        uncorrelated_variants_ch = GetUncorrelatedVariants(reference_bcf_files_ch)
            .collectFile(name: 'merged.prune.in', newLine: true).collect()

        // calculate the Z-scores for each parquet file
        zscore_ch = genes_ch
            .collate(100)
            .map{ gene_list -> CalculateZScore(permuted_parquet_ch, gene_list, uncorrelated_variants_ch) }
            .collectFile(name: 'pruned_z_scores.txt', skip: 1, keepHeader: true).collect()

        // Calculate gene gene matrix correlations
        uncorrelated_genes_ch = UncorrelatedGenes(zscore_ch, 0.1)

        // Get a collection of chunks for which to calculate LD
        loci_ch_raw = genes_ch
            .collate(100)
            .map{ gene_list -> CollectSignificantLoci(empirical_parquet_ch, gene_list)}
            .collectFile(name: 'loci_merged.txt', skip: 1, keepHeader: true).collect()

        // Add bp data to loci
        loci_annotated = AnnotateLoci(loci_ch_raw, variant_reference_ch, gene_reference_ch)

        // Flank loci and find the intersect between them
        loci_ch = IntersectLoci(
            loci_annotated.variant_loci, variant_flank_size,
            loci_annotated.gene_loci, gene_flank_size, genome_ref_ch)
            .splitText( by: 10 )

        // Calculate LD for all loci
        ld_ch = CalculateLd(permuted_parquet_ch, uncorrelated_genes_ch, loci_ch)

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
