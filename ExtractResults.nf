#!/usr/bin/env nextflow

/*
 * enables modules
 */
nextflow.enable.dsl = 2

// import modules
include { ExtractVariants; ExtractLociBed; ExtractLociAll } from './modules/CollectResults'
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
params.inclusion_step_output = 'NO_FILE'

if (params.help){
    helpmessage()
    exit 0
}

//Default parameters
Channel.fromPath(params.input).collect().set { input_parquet_ch }
Channel.fromPath(params.genes).splitCsv(header: ['gene']).map { row -> "${row.gene}" } .set { genes_ch }
Channel.fromPath(params.variant_reference).collect().set { variant_reference_ch }
Channel.fromPath(params.gene_reference).collect().set { gene_reference_ch }


inclusion_step_output_ch = file(params.inclusion_step_output)
bed_file_ch = file(params.bed)
variants_ch = file(params.variants)

Channel.fromPath(params.maf_table).set { maf_table_ch }

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

workflow {
    // Buffer genes
    genes_buffered_ch = genes_ch.collate(gene_chunk_size)

    if ( extract_variants ) {
        // Extract variants
        variants_extracted_ch = ExtractVariants(input_parquet_ch, variant_reference_ch, genes_buffered_ch, variants_ch, '+p_value,+z_score')
            .flatten()
            .map { file ->
                   def key = variants_ch.baseName
                   return tuple(key, file) }
            groupTuple()

        // Annotate loci
        variants_annotated_ch = AnnotateLoci(loci_extracted_ch, variant_reference_ch, gene_reference_ch, maf_table_ch, inclusion_dir_ch)

    }

    if ( extract_loci ) {
        // Extract loci
        loci_extracted_ch = ExtractLociAll(input_parquet_ch, loci_ch, variant_reference_ch, genes_buffered_ch, '+p_value,+z_score')
            .flatten()
            .map { file ->
                   def key = file.name.toString().tokenize('.').get(1)
                   return tuple(key, file) }
            groupTuple()

        // Annotate loci
        loci_annotated_ch = AnnotateLoci(loci_extracted_ch, variant_reference_ch, gene_reference_ch, maf_table_ch, inclusion_dir_ch)

    }

}

workflow.onComplete {
    println ( workflow.success ? "Pipeline finished!" : "Something crashed...debug!" )
}
