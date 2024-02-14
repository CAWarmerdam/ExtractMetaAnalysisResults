#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { ExtractSignificantResults; DefineFineMappingLoci } from './modules/CollectSignificantLoci'
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

params.maf_table = 'NO_FILE'
params.background_bed = 'NO_FILE'
params.inclusion_step_output = 'NO_FILE'

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

cohorts_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.cohort_new_name ]}
    .collect()

inclusion_step_output_ch = file(params.inclusion_step_output)
bed_file_ch = file(params.background_bed)

Channel.fromPath(params.maf_table).collect().set { maf_table_ch }

variant_flank_size=250000
gene_flank_size=1000000

gene_chunk_size=200
locus_chunk_size=100

enable_ld_calculation = true
enable_extract_loci = true
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

workflow GENE_CORRELATIONS {
    take:
        reference_bcf_files_ch
        permuted_parquet_ch
        variant_reference_ch
        genes_buffered_ch

    main:
        // Obtain a list of uncorrelated variants
        uncorrelated_variants_ch = GetUncorrelatedVariants(reference_bcf_files_ch)
            .collectFile(name: 'merged.prune.in', newLine: true, cache: 'lenient').collect()

        // Calculate the Z-scores for each gene list in the genes channel
        z_scores_split_ch = CalculateZScores(permuted_parquet_ch, variant_reference_ch, genes_buffered_ch, uncorrelated_variants_ch)

        // Combine Z-scores channel into a single file
        zscore_ch = z_scores_split_ch.collectFile(name: 'pruned_z_scores.txt', skip: 1, keepHeader: true, cache: 'lenient').collect()

        // Calculate gene gene matrix correlations
        uncorrelated_genes_out = UncorrelatedGenes(zscore_ch, 0.2)

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
        inclusion_step_output_ch
        genome_ref_ch
        cohorts_list

    main:
        // For every gene, get lead variants for significant results, and apply a window of 1Mb around the lead variant,
        // annotate with cis/trans write bed file
        significant_results_ch = ExtractSignificantResults(empirical_parquet_ch, genes_buffered_ch, variant_reference_ch, gene_reference_ch, inclusion_step_output_ch, 0.000000000002496, cohorts_list)
            .collectFile(name: 'lead_variants.csv', skip: 1, keepHeader: true, cache: true, storeDir: "${params.output}/significant_results").collect()

        // Output a channel of sorted bed files (by start pos).
        // Each bed file has the following columns: chr, start, end, name.
        // Where, each row is defined by a significant lead variant, with a 1Mb base pair window around the lead variant
        // The highest end pos must never be greater than the lowest start pos + 5Mb
        loci_ch = DefineFineMappingLoci(significant_results_ch, genome_ref_ch).flatten()
    emit:
        loci = loci_ch
}

workflow FINEMAPPING {
    take:
        empirical_parquet_ch
        permuted_parquet_ch
        variant_reference_ch
        uncorrelated_genes_ch
        loci_bed_ch

    main:
        ld_ch = CalculateLdMatrix(
                    permuted_parquet_ch, uncorrelated_genes_ch, variant_reference_ch,
                    loci_bed_ch)

        // Start finemapping here
}

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    // By default, always calculate gene correlations, and always run getting loci
    GENE_CORRELATIONS(reference_bcf_files_ch,permuted_parquet_ch,variant_reference_ch,genes_buffered_ch)

    uncorrelated_genes_buffered_ch = GENE_CORRELATIONS.out.uncorrelated_genes
        .splitCsv(header: ['gene']).map { row -> "${row.gene}" }.collate(gene_chunk_size)

    // ^^^ werkt
    // hieronder shaky

    // Define loci to do finemapping for
    LOCI(
        empirical_parquet_ch,genes_buffered_ch,
        variant_reference_ch,gene_reference_ch,genome_ref_ch,
        inclusion_step_output_ch,cohorts_ch.collect())

    // Do finemapping
    FINEMAPPING(
        empirical_parquet_ch,permuted_parquet_ch,variant_reference_ch,
        GENE_CORRELATIONS.out.uncorrelated_genes,LOCI.out.loci)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
