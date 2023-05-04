#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { ExtractSignificantResults; AnnotateResults; IntersectLoci; SelectFollowUpLoci; ExtractLoci; AnnotateLoci } from './modules/CollectSignificantLoci'
include { GetBreakpoints; CalculateAccuratePermutationPValues } from './modules/CalculateAccuratePermutationPValues'
include { CalculateLdMatrix; UncorrelatedGenes } from './modules/CalculateLdMatrix'
include { GetUncorrelatedVariants } from './modules/UncorrelatedVariants'
include { CalculateZScores } from './modules/CalculateZScores'
include { SampleOverlapMatrix } from './modules/Colocalization'


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

params.maf_table = 'NO_FILE'
params.clustered_loci = 'NO_FILE'
params.background_bed = 'NO_FILE'

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.empirical).collect().set { empirical_parquet_ch }
Channel.fromPath(params.permuted).collect().set { permuted_parquet_ch }
Channel.fromPath(params.reference_data).set { reference_bcf_files_ch }
Channel.fromPath(params.genes).splitCsv(header: ['gene']).map { row -> "${row.gene}" } .set { genes_ch }
Channel.fromPath(params.genome_reference).collect().set { genome_ref_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.gene_reference).collect().set { gene_reference_ch }
Channel.fromPath(params.inclusion_output).collect().set { inclusion_output_ch }

bed_file = file(params.background_bed)

Channel.fromPath(params.maf_table).set { maf_table_ch }

variant_flank_size=250000
gene_flank_size=250000

gene_chunk_size=20
locus_chunks_size=10

enable_ld_calculation = true
enable_extract_loci = false
enable_cis_trans_coloc = false

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


workflow ACCURATE_P_VALUES {
    // Obtain breakpoints to use for splitting variants
    breakpoints = GetBreakpoints(maf_table_ch.collect()).splitText()

    // For each of the breakpoints, get the uncorrelated variants
    uncorrelatedVariants = GetUncorrelatedVariants(breakpoints)
    CalculateAccuratePermutationPValues(
        empiricalResults, permutedResults,
        mafTable, breakpoints, uncorrelatedVariants, breakPoints, uncorrelatedVariants, andersonDarlingTable)

    CalculateAccuratePermutationPValues.out
}

workflow GENE_CORRELATIONS {
    take:
        reference_bcf_files_ch
        permuted_parquet_ch,
        variant_reference_ch,
        genes_buffered_ch

    main:
        // Obtain a list of uncorrelated variants
        uncorrelated_variants_ch = GetUncorrelatedVariants(reference_bcf_files_ch)
            .collectFile(name: 'merged.prune.in', newLine: true).collect()

        // Calculate the Z-scores for each gene list in the genes channel
        z_scores_split_ch = CalculateZScores(permuted_parquet_ch, variant_reference_ch, genes_buffered_ch, uncorrelated_variants_ch)

        // Combine Z-scores channel into a single file
        zscore_ch = z_scores_split_ch.collectFile(name: 'pruned_z_scores.txt', skip: 1, keepHeader: true).collect()

        // Calculate gene gene matrix correlations
        uncorrelated_genes_out = UncorrelatedGenes(zscore_ch, 0.05)

    emit:
        gene_correlations = uncorrelated_genes_out.correlations
        uncorrelated_genes = uncorrelated_genes_out.genes
}

workflow LOCI {
    take:
        empirical_parquet_ch
        genes_buffered_ch
        variant_reference_ch
        gene_reference_ch
        genome_ref_ch
        variant_flank_size
        gene_flank_size

    main:
        // Get a collection of chunks for which to calculate LD
        significant_results_ch = ExtractSignificantResults(empirical_parquet_ch, genes_buffered_ch, 0.000000000002496)
            .collectFile(name: 'loci_merged.txt', skip: 1, keepHeader: true).collect()

        // Add bp data to loci
        loci_bed_files = AnnotateResults(significant_results_ch, variant_reference_ch, gene_reference_ch)

        // Flank loci and find the union between them
        loci_ch = IntersectLoci(
            loci_bed_files.variant_loci, variant_flank_size,
            loci_bed_files.gene_loci, gene_flank_size, bed_file, genome_ref_ch)
    emit:
        loci = loci_ch
        loci_bed_files = loci_bed_files
}

workflow CALCULATE_LD {
    take:
        permuted_parquet_ch
        uncorrelated_genes_ch
        variant_reference_ch
        loci_ch
        locus_chunk_size

    main:
        // Calculate LD for all loci
        ld_ch = CalculateLdMatrix(
            permuted_parquet_ch, uncorrelated_genes_out.genes, variant_reference_ch,
            loci_ch.splitText( by: locus_chunk_size ))

    emit:
        ld_ch
}

workflow COLLECT_LOCI {
    take:
        empirical_parquet_ch
        genes_buffered_ch
        maf_table_ch
        loci_bed_files
        variant_flank_size
        gene_flank_size
        bed_file_ch
        genome_ref_ch
        locus_chunk_size

    main:

        // Flank loci and find the union between them, filtering on combinations where EITHER of the below is true:
        // 1. There is a cis genes involved
        // 2. There is a locus involved from the supplementary bed file
        loci_ch_pruned = FollowUpLoci(
            loci_bed_files.variant_loci, variant_flank_size,
            loci_bed_files.gene_loci, gene_flank_size, bed_file_ch, genome_ref_ch)

        // Extract empirical results for all significant loci
        loci_extracted_ch = ExtractLoci(empirical_parquet_ch, loci_ch_pruned, genes_buffered_ch)

        // Annotate loci
        loci_annotated_ch = AnnotateLoci(loci_extracted_ch, variant_reference_ch, gene_reference_ch, maf_table_ch)

    emit:
        loci_annotated_ch
}

workflow CIS_TRANS_COLOCALIZATION {
    take:
        loci_annotated_ch
        posterior_threshold
        cs_threshold
        output_cs_pip

    main:

        // For each cluster of overlapping significant loci, determine if the genes colocalize
        hypr_coloc_results = RunHyprColoc(
            loci_annotated_ch, posterior_threshold, cs_threshold, output_cs_pip)

    emit:
        hypr_coloc_results
}

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    // By default, always calculate gene correlations, and always run getting loci
    GENE_CORRELATIONS(reference_bcf_files_ch,permuted_parquet_ch,variant_reference_ch,genes_buffered_ch)

    // Extract significant results from the empirical side, and get loci as bed files
    LOCI(empirical_parquet_ch,genes_buffered_ch,variant_reference_ch,gene_reference_ch,genome_ref_ch,variant_flank_size,gene_flank_size)

    // If enabled run the following workflows:
    if ( enable_ld_calculation ) {
        CALCULATE_LD(permuted_parquet_ch,GENE_CORRELATIONS.uncorrelated_genes,variant_reference_ch,LOCI.loci,locus_chunk_size)
    }

    if ( enable_extract_loci ) {
        COLLECT_LOCI( empirical_parquet_ch,maf_table_ch,LOCI.loci_bed_files,variant_flank_size,gene_flank_size,bed_file_ch,genome_ref_ch,locus_chunk_size )
    }

    if ( enable_cis_trans_coloc ) {
        CIS_TRANS_COLOCALIZATION(COLLECT_LOCI.out,params.posterior_threshold,params.cs_threshold,params.output_cs_pip)
    }
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
