#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { SplitGeneVariantPairs; ExtractVariants; ExtractGeneVariantPairs; ExtractLociBed; ExtractLociAll } from './modules/CollectResults'
include { AnnotateLoci } from './modules/CollectSignificantLoci'



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
params.bed = 'NO_FILE'
params.variants = 'NO_FILE'
params.gene_variant_pairs = 'NO_FILE'
params.inclusion_step_output = 'NO_FILE'
params.cols = '+p_value,+z_score'
params.p_threshold = 'NULL'
params.output

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.input).collect().set { input_parquet_ch }
Channel.fromPath(params.permuted).collect().set { permuted_parquet_ch }
Channel.fromPath(params.genes).splitCsv(header: true).map { row -> "${row.ID}" } .set { genes_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.gene_reference).collect().set { gene_reference_ch }

cohorts_ch = Channel.fromPath(params.mastertable)
    .ifEmpty { error "Cannot find master table from: ${params.mastertable}" }
    .splitCsv(header: true, sep: '\t', strip: true)
    .map{row -> [ row.cohort_new_name ]}
    .collect()

bed_file_ch = file(params.bed)
variants_ch = file(params.variants)
gene_variant_pairs_ch = file(params.gene_variant_pairs)

Channel.fromPath(params.maf_table).collect().set { maf_table_ch }

gene_chunk_size=100
locus_chunk_size=100

extract_variants = params.variants != "NO_FILE"
extract_gene_variant_pairs = params.gene_variant_pairs != "NO_FILE"
extract_loci = params.bed != "NO_FILE"
annotate = params.maf_table != "NO_FILE"

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
summary['Extract loci']                             = extract_loci
summary['Extract variants']                         = extract_variants

log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "======================================================="

// Procedure:
// Gene-gene correlations

workflow {

    // Calculate gene-gene matrix correlations
    global_components = GlobalPrincipalExpressionPca(permuted_unbiased_snps_ch, n_threshold, 1000)

    // We need to calculate component effects on variants in a chunked manner.
    // Here, define chunks. We run a single job for each chunk
    variants_buffered_ch = variants.collate(variant_chunk_size)

    // We need to correct eQTL effects per
    // Here, define chunks. We run a single job for each chunk
    genes_buffered_ch = genes.collate(gene_chunk_size)

    // Estimate component effects on variants
    betas_ch = EstimateComponentEffectsOnVariants(input_parquet_ch, global_components, gene_reference_ch, variant_reference_ch, variants_buffered_ch, uncorrelated_genes)

    // Estimate effects corrected for k components
    corrected_eqtls_ch = CorrectQtlEffects(input_parquet_ch, global_components, betas_ch, genes_buffered_ch, pruned_variants)
}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
