#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { UncorrelatedGenes } from './modules/RunFineMapping'
include { CalculateZScoreMatrix; ConcatMatrix as ConcatZScoresMatrix; ConcatMatrix as ConcatSampleSizeMatrix } from './modules/CalculateZScores'

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
cohorts_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map { row -> [row.cohort_new_name] }

n_threshold = Channel.fromPath(new File(params.inclusion_dir, "filter_logs.log").toString())
    .splitCsv(header: true, sep: "\t", strip: true)
    .map { row -> [row.Dataset, row.N] }
    .join(cohorts_ch).view().map { row -> row[1] }.toInteger().sum().view()
Channel.fromPath(params.permuted).collect().set { permuted_parquet_ch }
Channel.fromPath(params.genes).splitCsv(header: true).map { row -> "${row.ID}" } .set { genes_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.uncorrelated_variants).collect().set { uncorrelated_variants_ch }

gene_chunk_size = 100

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

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    // Retrieve Z scores
    z_scores_split_ch = CalculateZScoreMatrix(permuted_parquet_ch, variant_reference_ch, genes_buffered_ch, uncorrelated_variants_ch, n_threshold, cohorts_ch.collect())

    // Combine Z-scores channel into a single file
    zscore_ch = ConcatZScoresMatrix(z_scores_split_ch.z_scores.collect(), Channel.value('mat.z.txt'))
    sample_size_ch = ConcatSampleSizeMatrix(z_scores_split_ch.sample_size.collect(), Channel.value('mat.n.txt'))

    zscore_ch.view()

    // Calculate gene gene matrix correlations
    uncorrelated_genes_out = UncorrelatedGenes(zscore_ch, sample_size_ch, n_threshold, 0.2)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
